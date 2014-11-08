// Copyright 2008 Google Inc.
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

package com.google.android.stardroid.units;

import android.util.FloatMath;

import com.google.android.stardroid.provider.ephemeris.OrbitalElements;
import com.google.android.stardroid.util.Geometry;
import com.google.android.stardroid.util.MathUtil;

/**
 * This class corresponds to an object's location in Euclidean space
 * when it is projected onto a unit sphere (with the earthbound viewer at the
 * center).
 *
 * @author Brent Bryan
 */
public class GeocentricCoordinates extends Vector3 {

  public GeocentricCoordinates(float x, float y, float z) {
    super(x, y, z);
  }

  /**
   * Subtracts the values of the given geocentric coordinates from this
   * object.
   */
  public void Subtract(GeocentricCoordinates other) {
    this.x -= other.x;
    this.y -= other.y;
    this.z -= other.z;
  }

  /** Recomputes the x, y, and z variables in this class based on the specified
   * {@link RaDec}.
   */
  public void updateFromRaDec(RaDec raDec) {
    updateFromRaDec(raDec.ra, raDec.dec);
  }

  private void updateFromRaDec(float ra, float dec) {
    float raRadians = ra * Geometry.DEGREES_TO_RADIANS;
    float decRadians = dec * Geometry.DEGREES_TO_RADIANS;

    this.x = FloatMath.cos(raRadians) * FloatMath.cos(decRadians);
    this.y = FloatMath.sin(raRadians) * FloatMath.cos(decRadians);
    this.z = FloatMath.sin(decRadians);
  }

  /**
   * Convert ra and dec to x,y,z where the point is place on the unit sphere.
   */
  public static GeocentricCoordinates getInstance(RaDec raDec) {
    return getInstance(raDec.ra, raDec.dec);
  }

  public static GeocentricCoordinates getInstance(float ra, float dec) {
    GeocentricCoordinates coords = new GeocentricCoordinates(0.0f, 0.0f, 0.0f);
    coords.updateFromRaDec(ra, dec);
    return coords;
  }

  /**
   * Convert ra and dec to x,y,z where the point is place on the unit sphere.
   */
  public static GeocentricCoordinates getInstanceFromFloatArray(float[] xyz) {
    return new GeocentricCoordinates(xyz[0], xyz[1], xyz[2]);
  }

  @Override
  public float[] toFloatArray() {
    return new float[] {x, y, z};
  }

  /**
   * Assumes it's an array of length 3.
   * @param xyz
   */
  public void updateFromFloatArray(float[] xyz) {
    this.x = xyz[0];
    this.y = xyz[1];
    this.z = xyz[2];
  }

  public void updateFromVector3(Vector3 v) {
    this.x = v.x;
    this.y = v.y;
    this.z = v.z;
  }

  @Override
  public GeocentricCoordinates copy() {
    return new GeocentricCoordinates(x, y, z);
  }

  public static GeocentricCoordinates getInstanceFromVector3(Vector3 v) {
    return new GeocentricCoordinates(v.x, v.y, v.z);
  }

  public static GeocentricCoordinates getInstance(OrbitalElements elem,
      GeocentricCoordinates zenith) {
    float anomaly = elem.getAnomaly();
    float ecc = elem.eccentricity;
    float radius = elem.distance * (1 - ecc * ecc) / (1 + ecc * MathUtil.cos(anomaly));

    // geocentric rectangular coordinates of satellite
    float per = elem.perihelion;
    float asc = elem.ascendingNode;
    float inc = elem.inclination;
    float xh = radius *
        (MathUtil.cos(asc) * MathUtil.cos(anomaly + per - asc) -
         MathUtil.sin(asc) * MathUtil.sin(anomaly + per - asc) *
         MathUtil.cos(inc));
    float yh = radius *
        (MathUtil.sin(asc) * MathUtil.cos(anomaly + per - asc) +
        MathUtil.cos(asc) * MathUtil.sin(anomaly + per - asc) *
        MathUtil.cos(inc));
    float zh = radius * (MathUtil.sin(anomaly + per - asc) * MathUtil.sin(inc));

    double EARTH_RADIUS_M = 6.371e6;  // meters
		double METERS_IN_AU = 149597870700d;
    double EARTH_RADIUS_AU = EARTH_RADIUS_M / METERS_IN_AU;
    GeocentricCoordinates viewer_coords = zenith;
    viewer_coords.normalize();
    viewer_coords.scale((float)EARTH_RADIUS_AU);  // Ignores deviation due to mountains, oblateness, etc.

    GeocentricCoordinates object_coords = new GeocentricCoordinates(xh, yh, zh);
    object_coords.Subtract(viewer_coords);
    object_coords.normalize();
    return object_coords;
  }
}