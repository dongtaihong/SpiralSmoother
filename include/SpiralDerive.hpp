/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-21 18:37:19
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-19 16:46:32
 * @FilePath: /SpiralSmoother/include/SpiralDerive.hpp
 * @Description:这个文件用于螺旋线的曲线系数求解、各项系数的偏导、曲线的各阶导
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include <cmath>
#include <vector>
using namespace std;

vector<double> coef_(6, 0);
vector<vector<double>> coef_deriv_(6, vector<double>(7, 0));

/**
 * @description: 各项系数的表达式
 * @return {*}
 */
void CoefInitialize(double thetai, double dthetai, double ddthetai,
                    double thetai_1, double dthetai_1, double ddthetai_1,
                    double s) {
  double s2 = s * s;
  double s3 = s2 * s;
  double s4 = s3 * s;
  double s5 = s2 * s3;
  double x0 = thetai;
  double dx0 = dthetai;
  double ddx0 = ddthetai;
  double x1 = thetai_1;
  double dx1 = dthetai_1;
  double ddx1 = ddthetai_1;
  coef_[0] = -6.0 * x0 / s5 - 3.0 * dx0 / s4 - 0.5 * ddx0 / s3 + 6.0 * x1 / s5 -
             3.0 * dx1 / s4 + 0.5 * ddx1 / s3;
  coef_[1] = 15.0 * x0 / s4 + 8.0 * dx0 / s3 + 1.5 * ddx0 / s2 -
             15.0 * x1 / s4 + 7.0 * dx1 / s3 - ddx1 / s2;
  coef_[2] = -10.0 * x0 / s3 - 6.0 * dx0 / s2 - 1.5 * ddx0 / s +
             10.0 * x1 / s3 - 4.0 * dx1 / s2 + 0.5 * ddx1 / s;
  coef_[3] = 0.5 * ddx0;
  coef_[4] = dx0;
  coef_[5] = x0;
}

/**      // plt::clf();
      // std::map<std::string, std::string> keywords;
      // keywords.insert(pair<string, string>("color", "r"));
      // keywords.insert(pair<string, string>("label", "平滑后路径"));
      // plt::plot(X, Y, keywords);

      // std::map<std::string, std::string> keywords1;
      // keywords1.insert(pair<string, string>("label", "原始路径点"));
      // keywords1.insert(pair<string, string>("marker", "p"));
      // plt::scatter(X2, Y2, 10.0, keywords1);

      // std::map<std::string, std::string> keywords2;
      // keywords2.insert(pair<string, string>("label", "螺旋线端点"));
      // keywords2.insert(pair<string, string>("marker", "^"));
      // plt::scatter(opt_xi, opt_yi, 10.0, keywords2);

      // std::map<std::string, std::string> keywords3;
      // keywords3.insert(pair<string, string>("label", "Spiral Smoother"));
      // plt::legend(keywords3);
      // plt::pause(0.05);
 * @return {*}
 */
void CoefDeriveInitialize(double thetai, double dthetai, double ddthetai,
                          double thetai_1, double dthetai_1, double ddthetai_1,
                          double s) {
  double s2 = s * s;
  double s3 = s2 * s;
  double s4 = s3 * s;
  double s5 = s2 * s3;
  double s6 = s3 * s3;
  double x0 = thetai;
  double dx0 = dthetai;
  double ddx0 = ddthetai;
  double x1 = thetai_1;
  double dx1 = dthetai_1;
  double ddx1 = ddthetai_1;
  // derive a
  // double a = -6.0 * x0 / p5 - 3.0 * dx0 / p4 - 0.5 * ddx0 / p3 + 6.0 * x1 /
  // p5 - 3.0 * dx1 / p4 + 0.5 * ddx1 / p3;
  coef_deriv_[0][0] = -6.0 / s5;
  coef_deriv_[0][1] = -3.0 / s4;
  coef_deriv_[0][2] = -0.5 / s3;
  coef_deriv_[0][3] = 6.0 / s5;
  coef_deriv_[0][4] = -3.0 / s4;
  coef_deriv_[0][5] = 0.5 / s3;
  coef_deriv_[0][6] = 30.0 * x0 / s6 + 12.0 * dx0 / s5 + 1.5 * ddx0 / s4 -
                      30.0 * x1 / s6 + 12.0 * dx1 / s5 - 1.5 * ddx1 / s4;

  // derive b
  // double b = 15.0 * x0 / p4 + 8.0 * dx0 / p3 + 1.5 * ddx0 / p2 - 15.0 * x1 /
  // p4 + 7.0 * dx1 / p3 - ddx1 / p2;
  coef_deriv_[1][0] = 15.0 / s4;
  coef_deriv_[1][1] = 8.0 / s3;
  coef_deriv_[1][2] = 1.5 / s2;
  coef_deriv_[1][3] = -15.0 / s4;
  coef_deriv_[1][4] = 7.0 / s3;
  coef_deriv_[1][5] = -1.0 / s2;
  coef_deriv_[1][6] = -60.0 * x0 / s5 - 24.0 * dx0 / s4 - 3.0 * ddx0 / s3 +
                      60.0 * x1 / s5 - 21.0 * dx1 / s4 + 2.0 * ddx1 / s3;

  // derive c
  // double c = -10.0 * x0 / p3 - 6.0 * dx0 / p2 - 1.5 * ddx0 / p + 10.0 * x1 /
  // p3 - 4.0 * dx1 / p2 + 0.5 * ddx1 / p;
  coef_deriv_[2][0] = -10.0 / s3;
  coef_deriv_[2][1] = -6.0 / s2;
  coef_deriv_[2][2] = -1.5 / s;
  coef_deriv_[2][3] = 10.0 / s3;
  coef_deriv_[2][4] = -4.0 / s2;
  coef_deriv_[2][5] = 0.5 / s;
  coef_deriv_[2][6] = 30.0 * x0 / s4 + 12.0 * dx0 / s3 + 1.5 * ddx0 / s2 -
                      30.0 * x1 / s4 + 8.0 * dx1 / s3 - 0.5 * ddx1 / s2;

  // derive d
  // double d = 0.5 * ddx0;
  coef_deriv_[3][2] = 0.5;

  // derive e
  // double e = dx0;
  coef_deriv_[4][1] = 1.0;

  // derive f
  // double f = x0;
  coef_deriv_[5][0] = 1.0;
}

/**
 * @description: 曲线各阶导的值
 * @param {int} type
 * @param {double} ratio
 * @param {double} delta_s
 * @return {*}
 */
double SpiralCurve(int type, double ratio, double delta_s) {
  double s = ratio * delta_s;
  double s2 = s * s;
  double s3 = s2 * s;
  double s4 = s3 * s;
  double s5 = s2 * s3;
  if (type == 0) {
    double thetas = coef_[0] * s5 + coef_[1] * s4 + coef_[2] * s3 +
                    coef_[3] * s2 + coef_[4] * s + coef_[5];
    return thetas;
  }
  if (type == 1) {
    double dthetas = 5 * coef_[0] * s4 + 4 * coef_[1] * s3 + 3 * coef_[2] * s2 +
                     2 * coef_[3] * s + coef_[4];
    return dthetas;
  }
  if (type == 2) {
    double ddthetas = 20 * coef_[0] * s3 + 12 * coef_[1] * s2 +
                      6 * coef_[2] * s + 2 * coef_[3];
    return ddthetas;
  }
  return 0;
}

/**
 * @description: dtheta对于thetai, dthetai, ddthetai, thetai_1, dthetai_1,
 * ddthetai_1, delta_s的偏导
 * @param {int} type
 * @param {double} ratio
 * @param {double} delta_s
 * @return {*}
 */
double DThetaDerivative(int type, double ratio, double delta_s) {
  double s = delta_s * ratio;
  double s2 = s * s;
  double s3 = s2 * s;
  double s4 = s2 * s2;

  double derivative = 5.0 * coef_deriv_[0][type] * s4 +
                      4.0 * coef_deriv_[1][type] * s3 +
                      3.0 * coef_deriv_[2][type] * s2 +
                      2.0 * coef_deriv_[3][type] * s + coef_deriv_[4][type];

  if (type == 6) {
    derivative += 20.0 * coef_[5] * s3 * ratio + 12.0 * coef_[4] * s2 * ratio +
                  6.0 * coef_[3] * s * ratio + 2.0 * coef_[2] * ratio;
  }
  return derivative;
}

/**
 * @description: ddtheta对于thetai, dthetai, ddthetai, thetai_1, dthetai_1,
 * ddthetai_1, delta_s的偏导
 * @param {int} type
 * @param {double} ratio
 * @param {double} delta_s
 * @return {*}
 */
double DDThetaDerivative(int type, double ratio, double delta_s) {
  double s = delta_s * ratio;
  double s2 = s * s;
  double s3 = s2 * s;

  double derivative =
      20.0 * coef_deriv_[0][type] * s3 + 12.0 * coef_deriv_[1][type] * s2 +
      6.0 * coef_deriv_[2][type] * s + 2.0 * coef_deriv_[3][type];

  if (type == 6) {
    derivative += 60.0 * coef_[0] * s2 * ratio + 24.0 * coef_[1] * s * ratio +
                  6.0 * coef_[2] * ratio;
  }
  return derivative;
}

/**
 * @description: x约束函数对于各个变量的偏导
 * @param {int} type
 * @param {int} m
 * @param {double} delta_si
 * @return {*}
 */
double GaussLegendreXDerivative(int type, int m, double delta_si) {
  vector<double> epsiloni = {-0.90617984594, -0.53846931011, 0, 0.53846931011,
                             0.90617984594};
  vector<double> wi = {0.236926885, 0.4786286705, 0.56888888889, 0.4786286705,
                       0.236926885};
  double gauss_legendre = 0;
  for (int i = 0; i < m; i++) {
    double segment_s = 0.5 * delta_si * epsiloni[i] + 0.5 * delta_si;
    double theta_segment_s = SpiralCurve(0, 1, segment_s);
    gauss_legendre += -wi[i] * sin(theta_segment_s) *
                      (coef_deriv_[0][type] * pow(segment_s, 5) +
                       coef_deriv_[1][type] * pow(segment_s, 4) +
                       coef_deriv_[2][type] * pow(segment_s, 3) +
                       coef_deriv_[3][type] * pow(segment_s, 2) +
                       coef_deriv_[4][type] * segment_s + coef_deriv_[5][type]);
  }
  double one = 0;
  if (type == 6) {
    for (int i = 0; i < m; i++) {
      double segment_s = 0.5 * delta_si * epsiloni[i] + 0.5 * delta_si;
      double theta_segment_s = SpiralCurve(0, 1, segment_s);
      double ratio = 0.5 * (epsiloni[i] + 1);
      one += wi[i] * cos(theta_segment_s);
      gauss_legendre += -wi[i] * sin(theta_segment_s) *
                        (5 * coef_[0] * pow(segment_s, 4) * ratio +
                         4 * coef_[1] * pow(segment_s, 3) * ratio +
                         3 * coef_[2] * pow(segment_s, 2) * ratio +
                         2 * coef_[3] * segment_s * ratio + coef_[4] * ratio);
    }
    one *= 0.5;
  }
  double derivative = one + 0.5 * delta_si * gauss_legendre;
  return derivative;
}

/**
 * @description: y约束函数对于各个变量的偏导
 * @param {int} type
 * @param {int} m
 * @param {double} delta_si
 * @return {*}
 */
double GaussLegendreYDerivative(int type, int m, double delta_si) {
  vector<double> epsiloni = {-0.90617984594, -0.53846931011, 0, 0.53846931011,
                             0.90617984594};
  vector<double> wi = {0.236926885, 0.4786286705, 0.56888888889, 0.4786286705,
                       0.236926885};
  double gauss_legendre = 0;
  for (int i = 0; i < m; i++) {
    double segment_s = 0.5 * delta_si * epsiloni[i] + 0.5 * delta_si;
    double theta_segment_s = SpiralCurve(0, 1, segment_s);
    gauss_legendre += wi[i] * cos(theta_segment_s) *
                      (coef_deriv_[0][type] * pow(segment_s, 5) +
                       coef_deriv_[1][type] * pow(segment_s, 4) +
                       coef_deriv_[2][type] * pow(segment_s, 3) +
                       coef_deriv_[3][type] * pow(segment_s, 2) +
                       coef_deriv_[4][type] * segment_s + coef_deriv_[5][type]);
  }
  double one = 0;
  if (type == 6) {
    for (int i = 0; i < m; i++) {
      double segment_s = 0.5 * delta_si * epsiloni[i] + 0.5 * delta_si;
      double theta_segment_s = SpiralCurve(0, 1, segment_s);
      double ratio = 0.5 * (epsiloni[i] + 1);
      one += wi[i] * sin(theta_segment_s);
      gauss_legendre += wi[i] * cos(theta_segment_s) *
                        (5 * coef_[0] * pow(segment_s, 4) * ratio +
                         4 * coef_[1] * pow(segment_s, 3) * ratio +
                         3 * coef_[2] * pow(segment_s, 2) * ratio +
                         2 * coef_[3] * segment_s * ratio + coef_[4] * ratio);
    }
    one *= 0.5;
  }
  double derivative = one + 0.5 * delta_si * gauss_legendre;
  return derivative;
}