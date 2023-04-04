/**
 * @brief Quaternion to Euler angles conversion
 * @cite Bernardes, Viollet - Quaternion to Euler angles conversion:
 *       A direct, general and computationally efficient method (2022)
 * @author Rishav
 * @date 2023-04-04
 */

/*
  compile: gcc main.c -o main -lm
  run: ./main
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846   // pi
#define PI_2 1.57079632679489661923 // pi/2

enum EulerSequence
{
  ZYZ = 323,
  ZXZ = 313,
  XYX = 121,
  XZX = 131,
  YXY = 212,
  YZY = 232,
  ZYX = 321,
  ZXY = 312,
  XYZ = 123,
  XZY = 132,
  YXZ = 213,
  YZX = 231
};

void parse_euler_sequence(const enum EulerSequence es, unsigned int *ijk)
{
  unsigned int sq = (int)es;
  ijk[0] = (sq / 100) % 10;
  ijk[1] = (sq / 10) % 10;
  ijk[2] = sq % 10;
}

// q = q0 + q1 * i + qz * j + q3 * k
// e = e1, e2, e3
double quat_to_euler(const double *q, const enum EulerSequence es, double *e)
{
  double tolerance = 1e-5;

  int n[3] = {0}; // ijk sequence
  unsigned int proper_flag = 1;
  parse_euler_sequence(es, n);

  if (n[0] == n[2])
  {
    n[2] = 6 - n[0] - n[1];
    proper_flag = 0;
  }

  double epsilon = (n[0] - n[1]) * (n[1] - n[2]) * (n[2] - n[1]) / 2.0f;

  double a, b, c, d;
  if (!proper_flag)
  {
    a = q[0] - q[n[1]];
    b = q[n[0]] + q[n[2]] * epsilon;
    c = q[n[1]] + q[0];
    d = q[n[2]] * epsilon - q[n[0]];
  }
  else
  {
    a = q[0];
    b = q[1];
    c = q[2];
    d = q[3] * epsilon;
  }

  e[1] = acos(2 * ((a * a + b * b) / (a * a + b * b + c * c + d * d)) - 1);
  double theta_plus = atan2(b, a);
  double theta_minus = atan2(d, c);

  if (abs(e[1]) < tolerance)
  {
    e[0] = 0.0f;
    e[2] = 2 * theta_plus - e[0];
  }
  else if (abs(e[1] - PI_2) < tolerance)
  {
    e[0] = 0.0f;
    e[2] = 2 * theta_plus + e[0];
  }
  else
  {
    e[0] = theta_plus - theta_minus;
    e[2] = theta_plus  + theta_minus;
  }

  if (!proper_flag)
  {
    e[2] = epsilon * e[2];
    e[1] = e[1] - PI_2;
  }
}

int main()
{
  double q[4] = {0.0, 0.0, 1.0, 0.0};
  double e[3] = {0.0f};

  quat_to_euler(q, ZYX, e);
  printf("q: [%f, %f, %f, %f]\r\n", q[0], q[1], q[2], q[3]);
  printf("e1: %f, e2: %f, e3: %f\r\n", e[0], e[1], e[2]);

  return 0;
}
