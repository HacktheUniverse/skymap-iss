// Copyright 2010 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package com.google.android.stardroid.layers;

import com.google.android.stardroid.R;
import com.google.android.stardroid.base.Closeables;
import com.google.android.stardroid.base.Lists;
import com.google.android.stardroid.base.TimeConstants;
import com.google.android.stardroid.control.AstronomerModel;
import com.google.android.stardroid.provider.ephemeris.OrbitalElements;
import com.google.android.stardroid.renderer.RendererObjectManager.UpdateType;
import com.google.android.stardroid.source.AbstractAstronomicalSource;
import com.google.android.stardroid.source.AstronomicalSource;
import com.google.android.stardroid.source.PointSource;
import com.google.android.stardroid.source.Sources;
import com.google.android.stardroid.source.TextSource;
import com.google.android.stardroid.source.impl.PointSourceImpl;
import com.google.android.stardroid.source.impl.TextSourceImpl;
import com.google.android.stardroid.units.GeocentricCoordinates;
import com.google.android.stardroid.util.Blog;
import com.google.android.stardroid.util.Geometry;
import com.google.android.stardroid.util.MiscUtil;
import com.google.android.stardroid.util.TimeUtil;

import android.content.res.Resources;
import android.graphics.Color;
import android.util.Log;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Date;
import java.util.EnumSet;
import java.util.List;
import java.util.TimeZone;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.TimeUnit;

/**
 * @author Brent Bryan
 */
public class IssLayer extends AbstractSourceLayer {
  private final ScheduledExecutorService scheduler = Executors.newScheduledThreadPool(1);
  private final AstronomerModel model;

  private IssSource issSource;

  public IssLayer(Resources resources, AstronomerModel model) {
    super(resources, true);
    this.model = model;
  }

  @Override
  protected void initializeAstroSources(ArrayList<AstronomicalSource> sources) {
    this.issSource = new IssSource(model, getResources());
    sources.add(issSource);

    scheduler.scheduleAtFixedRate(
        new OrbitalElementsGrabber(issSource), 0, 60, TimeUnit.SECONDS);
  }

  @Override
  public int getLayerId() {
    return -110;
  }

  @Override
  protected int getLayerNameId() {
    // TODO(brent): Update to different preference
    return R.string.show_planets_pref;
  }

  /** Thread Runnable which parses the orbital elements out of the Url. */
  static class OrbitalElementsGrabber implements Runnable {
    private static final long UPDATE_FREQ_MS = TimeConstants.MILLISECONDS_PER_HOUR;
    private static final String TAG = MiscUtil.getTag(OrbitalElementsGrabber.class);
    private static final String URL_STRING = "http://spaceflight.nasa.gov/realdata/" +
        "sightings/SSapplications/Post/JavaSSOP/orbit/ISS/SVPOST.html";

    private final IssSource source;
    private long lastSuccessfulUpdateMs = -1L;

    public OrbitalElementsGrabber(IssSource source) {
      this.source = source;
    }

    Date parseVectorTime(BufferedReader in) throws IOException,ParseException {
      String s;
      while ((s = in.readLine()) != null && !s.contains("Vector Time (GMT)")) {}
      String date_string = s.trim().split("\\s+")[3];
      DateFormat df = new SimpleDateFormat("yyyy/DDD/kk:mm:ss.SSS");
      df.setTimeZone(TimeZone.getTimeZone("GMT"));
      return df.parse(date_string);
    }

    /**
     * Parses the OrbitalElements from the given BufferedReader.  Factored out
     * of {@link #getOrbitalElements} to simplify testing.
     */
    OrbitalElements parseOrbitalElements(BufferedReader in)
        throws IOException,ParseException {
      Date epoch = parseVectorTime(in);
      String s;
      while ((s = in.readLine()) != null && !s.contains("M50 Keplerian")) {}

      // Skip the dashed line
      in.readLine();

      float[] params = new float[9];
      int i = 0;
      for (; i < params.length && (s = in.readLine()) != null; i++) {
        s = s.substring(46).trim();
        String[] tokens = s.split("\\s+");
        params[i] = Float.parseFloat(tokens[2]);
      }


      if (i == params.length) {  // we read all the data.
        StringBuilder sb = new StringBuilder();
        for (int pi = 0; pi < params.length; pi++) {
          sb.append(" " + params[pi]);
        }
        Blog.d(this, "Params: " + sb);
        OrbitalElements elements = convertNasaDataToOrbitalElements(params);
        elements.epoch = epoch;
        return elements;
      }
      return null;
    }

    OrbitalElements convertNasaDataToOrbitalElements(float[] params) {
      // See http://spaceflight1.nasa.gov/realdata/elements/index.html
      float a_m = params[0];  // Semimajor axis in meters
      float e = params[1];  // Eccentricity
      float i = params[2];  // Inclination in degrees
      float Wp = params[3];  // lower-case-omega, argument of perigee in degrees
      // Right Ascension of the ascending node in degrees, a.k.a. upper-case-omega, a.k.a. Longitude of the ascending node.
      float RA = params[4];
      float TA = params[5];  // True anomaly in degrees
      float MA = params[6];  // Mean anomaly in degrees
      float Ha = params[7];  // height of apogee nautical miles
      float Hp = params[8];  // height of perigee, nautical miles

      // http://en.wikipedia.org/wiki/Longitude_of_the_periapsis#Calculation_from_state_vectors
      float longitude_of_perigee_degrees = Wp + RA;

      // http://en.wikipedia.org/wiki/Mean_longitude#Calculation
      float mean_longitude_degrees = MA + longitude_of_perigee_degrees;

      float METERS_IN_AU = 149597870700f;
      float semimajor_axis_au = a_m / METERS_IN_AU;
      float inclination_radians = i * Geometry.DEGREES_TO_RADIANS;
      float right_ascension_node_radians = RA * Geometry.DEGREES_TO_RADIANS;
      float longitude_of_perigee_radians = Wp * Geometry.DEGREES_TO_RADIANS;
      float mean_longitude_radians = mean_longitude_degrees * Geometry.DEGREES_TO_RADIANS;

      // NOTE: These are the M50 orbital elements, i.e. as of 1950-Jan-1.
      // TODO: Compensate for earth's precession between M50 (a.k.a. B1950) and
      // J2000 (which is used for all other calculations here).
      // See http://earth-info.nga.mil/GandG/publications/tr8350.2/tr8350.2-a/Appendix.pdf
      return new OrbitalElements(
          semimajor_axis_au,
          e,
          inclination_radians,
          right_ascension_node_radians,
          longitude_of_perigee_radians,
          mean_longitude_radians);
    }

    /**
     * Reads the given URL and returns the OrbitalElements associated with the object
     * described therein. NOTE: These are M50 geocentric orbital elements.
     */
    OrbitalElements getOrbitalElements(String urlString) {
      BufferedReader in = null;
      try {
        URLConnection connection = new URL(urlString).openConnection();
        in = new BufferedReader(new InputStreamReader(connection.getInputStream()));
        return parseOrbitalElements(in);
      } catch (IOException e) {
        Log.e(TAG, "Error reading Orbital Elements");
      } catch (ParseException e) {
        Log.e(TAG, "Error parsing Orbital Elements epoch");
      } finally {
        Closeables.closeSilently(in);
      }
      return null;
    }

    @Override
    public void run() {
      long currentTimeMs = System.currentTimeMillis();
      if ((currentTimeMs - lastSuccessfulUpdateMs) > UPDATE_FREQ_MS) {
        Blog.d(this, "Fetching ISS data...");
        OrbitalElements elements = getOrbitalElements(URL_STRING);
        if (elements == null) {
          Log.d(TAG, "Error downloading ISS orbital data");
        } else {
          lastSuccessfulUpdateMs = currentTimeMs;
          source.setOrbitalElements(elements);
        }
      }
    }
  }

  /** AstronomicalSource corresponding to the International Space Station. */
  static class IssSource extends AbstractAstronomicalSource {
    private static final long UPDATE_FREQ_MS = 1L * TimeConstants.MILLISECONDS_PER_SECOND;
    private static final int ISS_COLOR = Color.YELLOW;

    // TODO: Don't set a default position, so that the ISS doesn't appear in
    // a fixed location if you don't have a net connection.
    private final GeocentricCoordinates coords = new GeocentricCoordinates(1f, 0f, 0f);
    private final ArrayList<PointSource> pointSources = new ArrayList<PointSource>();
    private final ArrayList<TextSource> textSources = new ArrayList<TextSource>();
    private final AstronomerModel model;
    private final String name;

    private OrbitalElements orbitalElements = null;
    private boolean orbitalElementsChanged;
    private long lastUpdateTimeMs = 0L;

    public IssSource(AstronomerModel model, Resources resources) {
      this.model = model;
      this.name = resources.getString(R.string.space_station);

      pointSources.add(new PointSourceImpl(coords, ISS_COLOR, 5));
      textSources.add(new TextSourceImpl(coords, name, ISS_COLOR));
    }

    public synchronized void setOrbitalElements(OrbitalElements elements) {
      // TODO: Only set this to true if the elements actually changed.
      this.orbitalElements = elements;
      orbitalElementsChanged = true;
    }

    @Override
    public List<String> getNames() {
      // TODO: Change the search function to query getNames() more than once,
      // so that ISS can be left out of the index until its location is known.
      return Lists.asList(name);
    }

    @Override
    public GeocentricCoordinates getSearchLocation() {
      return coords;
    }

    private void updateCoords(Date time) {
      lastUpdateTimeMs = time.getTime();
      orbitalElementsChanged = false;

      if (orbitalElements == null) {
        return;
      }

      OrbitalElements elements_now;
      if (orbitalElements.epoch == null) {
        // |orbitalElements| is kept up to date.
        elements_now = orbitalElements;
      } else {
        // Project |orbitalElements| forward to now.
        // values from NASA factsheet: http://nssdc.gsfc.nasa.gov/planetary/factsheet/fact_notes.html
        double EARTH_MASS = 5.9726e24;  // kg
        double G = 6.67384e-11;  // meters^3 kilograms^-1 seconds^-2
        double METERS_IN_AU = 149597870700d;
        double a = orbitalElements.distance * METERS_IN_AU;
        double mu = G * EARTH_MASS;
        double period = 2 * Math.PI * Math.sqrt(a * a * a / mu);  // seconds

        long time_delta_ms = lastUpdateTimeMs - orbitalElements.epoch.getTime();
        double time_delta = ((double)time_delta_ms) / 1000;  // seconds
        double orbits = time_delta / period;
        double phase_radians = 2 * Math.PI * (orbits % 1.0);
        double mean_longitude_now = (orbitalElements.meanLongitude + phase_radians) % (2 * Math.PI);

        elements_now = new OrbitalElements(
          orbitalElements.distance,
          orbitalElements.eccentricity,
          orbitalElements.inclination,
          orbitalElements.ascendingNode,
          orbitalElements.perihelion,
          (float)mean_longitude_now
        );
      }

      coords.assign(GeocentricCoordinates.getInstance(elements_now, model.getZenith()));
      Blog.d(this, "Changed ISS coords to " + coords.toString());
    }

    @Override
    public Sources initialize() {
      updateCoords(model.getTime());
      return this;
    }

    @Override
    public synchronized EnumSet<UpdateType> update() {
      EnumSet<UpdateType> updateTypes = EnumSet.noneOf(UpdateType.class);

      Date modelTime = model.getTime();
      if (orbitalElementsChanged ||
          Math.abs(modelTime.getTime() - lastUpdateTimeMs) > UPDATE_FREQ_MS) {

        updateCoords(modelTime);
        if (orbitalElements != null) {
          updateTypes.add(UpdateType.UpdatePositions);
        }
      }
      return updateTypes;
    }

    @Override
    public List<? extends TextSource> getLabels() {
      // TODO: Change the display update to query getLabels() more than once,
      // so that ISS can be left out of the view until its location is known.
      return textSources;
    }

    @Override
    public List<? extends PointSource> getPoints() {
      // TODO: Change the display update to query getPoints() more than once,
      // so that ISS can be left out of the view until its location is known.
      return pointSources;
    }
  }


}
