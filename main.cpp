#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

const double PI = 3.141592653589793;
typedef vector<double> Vector;

Vector
operator* (const double &scalar, const Vector &vec)
{
  Vector result;
  for (double value : vec)
    {
      result.push_back (scalar * value);
    }
  return result;
}
Vector
operator/ (const Vector &vec, const double &scalar)
{
  Vector result;
  for (double value : vec)
    {
      result.push_back (value / scalar);
    }
  return result;
}
Vector
operator+ (const Vector &v1, const Vector &v2)
{
  Vector result;
  if (v1.size () != v2.size ())
    {
      cout << "Размеры векторов не совпадают!!!";
    }
  else
    {
      for (int i = 0; i < v1.size (); i++)
        {
          result.push_back (v1[i] + v2[i]);
        }
    }
  return result;
}
Vector
operator- (const Vector &v1, const Vector &v2)
{
  Vector result;
  if (v1.size () != v2.size ())
    {
      cout << "Размеры векторов не совпадают!!!";
    }
  else
    {
      for (int i = 0; i < v1.size (); i++)
        {
          result.push_back (v1[i] - v2[i]);
        }
    }
  return result;
}
double
abs (const Vector &v)
{
  double result = 0;
  for (double value : v)
    {
      if (abs (value) > result)
        result = abs (value);
    }
  return result;
}

pair<Vector, Vector>
RungeKutta (double h, double t, Vector x, Vector (*f) (double, const Vector &))
{
  Vector k1 = h * f (t, x);
  Vector k2 = h * f (t + h / 2, x + k1 / 2);
  Vector k3 = h * f (t + h / 2, x + k2 / 2);
  Vector k4 = h * f (t + h, x + k3);
  return { x + (k1 + 2 * k2 + 2 * k3 + k4) / 6, 2 * (k1 - k2 - k3 + k4) / 3 };
}

class ShootingMethod
{
private:
  Vector (*f) (double, const Vector &);
  double t0;
  double t1;
  double tol_runge_kutta;
  double tol_deriv;
  double tol_newton;
  Vector x1;

public:
  ShootingMethod (Vector (*f) (double, const Vector &), double t0, double t1,
                  double tol_runge_kutta, double tol_newton, double tol_deriv,
                  Vector x1)
      : f (f), t0 (t0), t1 (t1), tol_runge_kutta (tol_runge_kutta),
        tol_newton (tol_newton), tol_deriv (tol_deriv), x1 (x1)
  {
  }
  Vector
  integrate (const Vector &x0)
  {
    double h = (t1 - t0) * 1e-4;
    double t = t0;
    Vector x = x0;
    while (t <= t1)
      {
        pair<Vector, Vector> x_temp_error = RungeKutta (h, t, x, f);

        x = x_temp_error.first;
        t = t + h;
      }

    return x;
  }

  Vector
  newtonMethod (int maxit, Vector x0 = { 0, 0 })
  {
    Vector x = x0;
    for (int i = 0; i < maxit; i++)
      {
        Vector F_value = ShootingMethod::F (x);
        Vector dFda = dF_da (x);
        Vector dFdb = dF_db (x);

        double detJ = dFda[0] * dFdb[1] - dFda[1] * dFdb[0];

        double deltaX = (dFdb[1] * F_value[0] - dFdb[0] * F_value[1]) / detJ;
        double deltaY = (-dFda[1] * F_value[0] + dFda[0] * F_value[1]) / detJ;

        if (std::abs (deltaX) < tol_newton && std::abs (deltaY) < tol_newton)
          {

            break;
          }
        x[0] -= deltaX;
        x[1] -= deltaY;
      }
    return x;
  }
  Vector
  F (const Vector &x)
  {
    Vector result = integrate ({ 0, x[1], x[0], x[1] });
    return { result[0] - x1[0], result[3] - x1[1] }; // двумерный
  }

  Vector
  dF_da (const Vector &x)
  {
    Vector result
        = (F ({ x[0] + tol_deriv, x[1] }) - F ({ x[0] - tol_deriv, x[1] }))
          / (2 * tol_deriv);
    return result;
  }

  Vector
  dF_db (const Vector &x)
  {
    Vector result
        = (F ({ x[0], x[1] + tol_deriv }) - F ({ x[0], x[1] - tol_deriv }))
          / (2 * tol_deriv);
    return result;
  }
};

Vector
f (double t, const Vector &x)
{
  double alpha = 10.5;
  Vector dxdt
      = { x[1], x[3] - x[0] / (1 + alpha * x[0] * x[0]),
          x[3] * (1 - alpha * x[0] * x[0])
              / ((1 + alpha * x[0] * x[0]) * (1 + alpha * x[0] * x[0])),
          -x[2] };
  return dxdt;
}

int
main ()
{
  double t0 = 0;
  double t1 = PI * 0.5;
  double tol_runge_kutta = 1e-14;
  double tol_newton = 1e-6;
  double tol_deriv = 1e-6;
  Vector x1 = { 1, 0 };
  ShootingMethod sm (f, t0, t1, tol_runge_kutta, tol_newton, tol_deriv, x1);
  double maxit = 100;
  Vector x = sm.newtonMethod (maxit);

  cout << std::setprecision (10);
  for (double value : x)
    {
      cout << value << " ";
    }
  cout << endl;

  return 0;
}