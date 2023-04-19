/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-20 22:09:26
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-19 17:06:30
 * @FilePath: /SpiralSmoother/src/ipopt_interface.cpp
 * @Description: 螺旋线问题的非线性ipopt建模
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include "ipopt_interface.hpp"

#include "SpiralDerive.hpp"
#include "matplotlibcpp.h"

double weight_curve_length_ = 100.0;
double weight_kappa_ = 0.1;
double weight_dkappa_ = 5.0;

namespace plt = matplotlibcpp;

// result
vector<double> opt_thetai;
vector<double> opt_dthetai;
vector<double> opt_ddthetai;
vector<double> opt_xi;
vector<double> opt_yi;
vector<double> opt_deltas;

void CsvSourceData::GetSourceData(const char* filename) {
  ifstream readCSV(filename, ios::in);
  string line;
  int input_index = 0;
  int line_number = 0;
  double distance_ = 0;
  double heading_ = 0;
  while (getline(readCSV, line)) {
    stringstream ss(line);
    string s;
    vector<string> temp_line;
    pair<double, double> temp_point;
    while (getline(ss, s, ',')) {
      temp_line.emplace_back(s);
    }
    temp_point.first = stod(temp_line[0]);
    temp_point.second = stod(temp_line[1]);
    //将起点放入
    if (line_number == 0) {
      input_points.emplace_back(temp_point);
      line_number++;
      continue;
    }

    distance_ =
        sqrt(pow(temp_point.first - input_points[input_index].first, 2) +
             pow(temp_point.second - input_points[input_index].second, 2));
    heading_ = atan2((temp_point.second - input_points.back().second),
                     (temp_point.first - input_points.back().first));
    //每隔2.5米左右采一个参考点（间隔可以自己去调，比如1米就把这里判断条件改成1左右）
    if (distance_ > 1.7 && distance_ < 3.3) {
      input_points.emplace_back(temp_point);
      this->distance.emplace_back(distance_);
      this->relative_theta_.emplace_back(heading_);
      input_index++;
    }
    line_number++;
  }
  this->number = input_points.size();
  relative_theta_.emplace_back(relative_theta_.back());
}

// constructor
SpiralSmootherNLP::SpiralSmootherNLP() {}

// destructor
SpiralSmootherNLP::~SpiralSmootherNLP() {}

// returns the size of the problem
bool SpiralSmootherNLP::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                                     Ipopt::Index& nnz_jac_g,
                                     Ipopt::Index& nnz_h_lag,
                                     Ipopt::TNLP::IndexStyleEnum& index_style) {
  //目标变量
  n = SourceData.number * 5 + SourceData.number - 1;
  //约束函数规模
  m = (SourceData.number - 1) * 2 + SourceData.number;
  //约束函数jacobian矩阵非零量
  nnz_jac_g = (SourceData.number - 1) * 2 * 9 + SourceData.number * 2;
  //目标hessian矩阵非零量
  nnz_h_lag = 0;
  //索引风格
  index_style = IndexStyleEnum::C_STYLE;
  return true;
}

// returns the variable bounds
bool SpiralSmootherNLP::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l,
                                        Ipopt::Number* x_u, Ipopt::Index m,
                                        Ipopt::Number* g_l,
                                        Ipopt::Number* g_u) {
  //每次循环设置的是第i个点的五个量：theta, kappa, dkappa, x, y
  for (int i = 0; i < SourceData.number; ++i) {
    int index = i * 5;
    double theta_lower = 0.0;
    double theta_upper = 0.0;
    double kappa_lower = 0.0;
    double kappa_upper = 0.0;
    double dkappa_lower = 0.0;
    double dkappa_upper = 0.0;
    double x_lower = 0.0;
    double x_upper = 0.0;
    double y_lower = 0.0;
    double y_upper = 0.0;
    if (i == 0) {
      //起点约束
      theta_lower = SourceData.relative_theta_[0];
      theta_upper = SourceData.relative_theta_[0];
      kappa_lower = -0.0;
      kappa_upper = 0.0;
      dkappa_lower = -0.01;
      dkappa_upper = 0.01;
      x_lower = SourceData.input_points[0].first;
      x_upper = SourceData.input_points[0].first;
      y_lower = SourceData.input_points[0].second;
      y_upper = SourceData.input_points[0].second;
    } else if (i == SourceData.number - 1) {
      //终点约束
      theta_lower = SourceData.relative_theta_.back();
      theta_upper = SourceData.relative_theta_.back();
      kappa_lower = -0.0;
      kappa_upper = 0.0;
      dkappa_lower = -0.01;
      dkappa_upper = 0.01;
      x_lower = SourceData.input_points.back().first;
      x_upper = SourceData.input_points.back().first;
      y_lower = SourceData.input_points.back().second;
      y_upper = SourceData.input_points.back().second;
    } else {
      //中间点的运动学与边界约束
      theta_lower = SourceData.relative_theta_[i] - M_PI * 0.5;
      theta_upper = SourceData.relative_theta_[i] + M_PI * 0.5;
      kappa_lower = -0.10;
      kappa_upper = 0.10;
      dkappa_lower = -0.01;
      dkappa_upper = 0.01;
      x_lower = SourceData.input_points[i].first - SourceData.boundary;
      x_upper = SourceData.input_points[i].first + SourceData.boundary;
      y_lower = SourceData.input_points[i].second - SourceData.boundary;
      y_upper = SourceData.input_points[i].second + SourceData.boundary;
    }
    //将上面的上下界填到x_l, x_u中
    // theta
    x_l[index] = theta_lower;
    x_u[index] = theta_upper;
    // kappa
    x_l[index + 1] = kappa_lower;
    x_u[index + 1] = kappa_upper;
    // dkappa
    x_l[index + 2] = dkappa_lower;
    x_u[index + 2] = dkappa_upper;
    // x
    x_l[index + 3] = x_lower;
    x_u[index + 3] = x_upper;
    // y
    x_l[index + 4] = y_lower;
    x_u[index + 4] = y_upper;
  }

  // delta_s的边界约束
  int variable_offset = SourceData.number * 5;
  for (int i = 0; i < SourceData.number - 1; ++i) {
    x_l[variable_offset + i] =
        SourceData.distance[i] - 2.0 * SourceData.boundary;
    x_u[variable_offset + i] = SourceData.distance[i] * M_PI * 0.5;
  }

  //段连接点的位置约束：等式约束，上下界都是0
  for (int i = 0; i < SourceData.number - 1; ++i) {
    // x
    g_l[i * 2] = 0.0;
    g_u[i * 2] = 0.0;
    // y
    g_l[i * 2 + 1] = 0.0;
    g_u[i * 2 + 1] = 0.0;
  }
  //每个点的位置约束
  int constraint_offset = 2 * (SourceData.number - 1);
  for (int i = 0; i < SourceData.number; ++i) {
    g_l[constraint_offset + i] = 0.0;
    g_u[constraint_offset + i] = pow(SourceData.boundary, 2);
  }
  return true;
}

// returns the initial point for the problem
bool SpiralSmootherNLP::get_starting_point(Ipopt::Index n, bool init_x,
                                           Ipopt::Number* x, bool init_z,
                                           Ipopt::Number* z_L,
                                           Ipopt::Number* z_U, Ipopt::Index m,
                                           bool init_lambda,
                                           Ipopt::Number* lambda) {
  //根据原始数据来赋值第i个点的theta, kappa, dkappa, x, y
  for (int i = 0; i < SourceData.number; ++i) {
    int index = i * 5;
    x[index] = SourceData.relative_theta_[i];
    x[index + 1] = 0.0;
    x[index + 2] = 0.0;
    x[index + 3] = SourceData.input_points[i].first;
    x[index + 4] = SourceData.input_points[i].second;
  }

  //根据原始数据来赋值第i段的delta_s
  int variable_offset = SourceData.number * 5;
  for (int i = 0; i < SourceData.number - 1; ++i) {
    double delta_theta =
        SourceData.relative_theta_[i + 1] - SourceData.relative_theta_[i];
    x[variable_offset + i] = SourceData.distance[i] / cos(0.5 * delta_theta);
  }

  //曲率的初始值
  for (int i = 0; i < SourceData.number - 1; ++i) {
    double delta_theta =
        SourceData.relative_theta_[i + 1] - SourceData.relative_theta_[i];
    x[(i + 1) * 5 + 1] = delta_theta / x[variable_offset + i];
  }
  x[1] = x[6];
  return true;
}

// returns the value of the objective function
bool SpiralSmootherNLP::eval_f(Ipopt::Index n, const Ipopt::Number* x,
                               bool new_x, Ipopt::Number& obj_value) {
  obj_value = 0.0;
  for (int i = 0; i < SourceData.number - 1; ++i) {
    double delta_si = x[5 * SourceData.number + i];
    double theta_i = x[5 * i + 0];
    double dtheta_i = x[5 * i + 1];
    double ddtheta_i = x[5 * i + 2];
    double theta_i_1 = x[5 * (i + 1) + 0];
    double dtheta_i_1 = x[5 * (i + 1) + 1];
    double ddtheta_i_1 = x[5 * (i + 1) + 2];
    CoefInitialize(theta_i, dtheta_i, ddtheta_i, theta_i_1, dtheta_i_1,
                   ddtheta_i_1, delta_si);
    //第一项
    obj_value += delta_si * weight_curve_length_;

    //从每段螺旋线中，采m个点，然后把每个点的kappa加起来，整体越平滑的线，kappa之和越小
    for (int j = 0; j < SourceData.m; ++j) {
      double ratio = static_cast<double>(j) / static_cast<double>(SourceData.m);

      //第二项
      double kappa = SpiralCurve(1, ratio, delta_si);
      obj_value += kappa * kappa * weight_kappa_;

      //第三项
      double dkappa = SpiralCurve(2, ratio, delta_si);
      obj_value += dkappa * dkappa * weight_dkappa_;
    }
  }
  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool SpiralSmootherNLP::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x,
                                    bool new_x, Ipopt::Number* grad_f) {
  int variable_offset = SourceData.number * 5;
  for (int i = 0; i + 1 < SourceData.number; ++i) {
    int index0 = i * 5;
    int index1 = (i + 1) * 5;
    double delta_si = x[5 * SourceData.number + i];

    grad_f[variable_offset + i] = weight_curve_length_ * 1.0;
    for (int j = 0; j < SourceData.m; ++j) {
      double ratio = static_cast<double>(j) / static_cast<double>(SourceData.m);

      double kappa = SpiralCurve(1, ratio, delta_si);
      grad_f[index0] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(0, ratio, delta_si);
      grad_f[index0 + 1] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(1, ratio, delta_si);
      grad_f[index0 + 2] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(2, ratio, delta_si);
      grad_f[index1] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(3, ratio, delta_si);
      grad_f[index1 + 1] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(4, ratio, delta_si);
      grad_f[index1 + 2] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(5, ratio, delta_si);
      grad_f[variable_offset + i] +=
          weight_kappa_ * 2.0 * kappa * DThetaDerivative(6, ratio, delta_si);

      double dkappa = SpiralCurve(2, ratio, delta_si);
      grad_f[index0] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(0, ratio, delta_si);
      grad_f[index0 + 1] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(1, ratio, delta_si);
      grad_f[index0 + 2] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(2, ratio, delta_si);
      grad_f[index1] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(3, ratio, delta_si);
      grad_f[index1 + 1] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(4, ratio, delta_si);
      grad_f[index1 + 2] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(5, ratio, delta_si);
      grad_f[variable_offset + i] +=
          weight_dkappa_ * 2.0 * dkappa * DDThetaDerivative(6, ratio, delta_si);
    }
  }
  return true;
}

// return the value of the constraints: g(x)
bool SpiralSmootherNLP::eval_g(Ipopt::Index n, const Ipopt::Number* x,
                               bool new_x, Ipopt::Index m, Ipopt::Number* g) {
  //公式1 2，分别N-1维
  for (int i = 0; i + 1 < SourceData.number; ++i) {
    int index0 = i * 5;
    int index1 = (i + 1) * 5;
    double delta_si = x[5 * SourceData.number + i];
    double theta_i = x[5 * i + 0];
    double dtheta_i = x[5 * i + 1];
    double ddtheta_i = x[5 * i + 2];
    double theta_i_1 = x[5 * (i + 1) + 0];
    double dtheta_i_1 = x[5 * (i + 1) + 1];
    double ddtheta_i_1 = x[5 * (i + 1) + 2];
    CoefInitialize(theta_i, dtheta_i, ddtheta_i, theta_i_1, dtheta_i_1,
                   ddtheta_i_1, delta_si);

    //这一部分是求gauss-legendre数值积分
    vector<double> epsiloni = {-0.90617984594, -0.53846931011, 0, 0.53846931011,
                               0.90617984594};
    vector<double> wi = {0.236926885, 0.4786286705, 0.56888888889, 0.4786286705,
                         0.236926885};
    double gauss_legendre_x = 0;
    double gauss_legendre_y = 0;
    for (int j = 0; j < SourceData.m; ++j) {
      double segment_s = 0.5 * delta_si * epsiloni[j] + 0.5 * delta_si;
      double theta_segment_s = SpiralCurve(0, 1, segment_s);
      gauss_legendre_x += wi[j] * cos(theta_segment_s);
      gauss_legendre_y += wi[j] * sin(theta_segment_s);
    }
    gauss_legendre_x *= delta_si / 2;
    gauss_legendre_y *= delta_si / 2;

    double x_diff = x[index1 + 3] - x[index0 + 3] - gauss_legendre_x;
    g[i * 2] = x_diff;
    double y_diff = x[index1 + 4] - x[index0 + 4] - gauss_legendre_y;
    g[i * 2 + 1] = y_diff;
  }
  //公式3，N维
  int constraint_offset = 2 * (SourceData.number - 1);
  for (int i = 0; i < SourceData.number; ++i) {
    int variable_index = i * 5;
    double x_cor = x[variable_index + 3];
    double y_cor = x[variable_index + 4];

    double x_diff = x_cor - SourceData.input_points[i].first;
    double y_diff = y_cor - SourceData.input_points[i].second;

    g[constraint_offset + i] = x_diff * x_diff + y_diff * y_diff;
  }
  return true;
}

// return the structure or values of the jacobian
bool SpiralSmootherNLP::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x,
                                   bool new_x, Ipopt::Index m,
                                   Ipopt::Index nele_jac, Ipopt::Index* iRow,
                                   Ipopt::Index* jCol, Ipopt::Number* values) {
  if (values == nullptr) {
    int nz_index = 0;

    int variable_offset = SourceData.number * 5;
    for (int i = 0; i + 1 < SourceData.number; ++i) {
      int variable_index = i * 5;

      // theta0
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 0;
      ++nz_index;

      // kappa0
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 1;
      ++nz_index;

      // dkappa0
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 2;
      ++nz_index;

      // x0
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 3;
      ++nz_index;

      // theta1
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 5;
      ++nz_index;

      // kappa1
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 6;
      ++nz_index;

      // dkappa1
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 7;
      ++nz_index;

      // x1
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_index + 8;
      ++nz_index;

      // s
      iRow[nz_index] = i * 2;
      jCol[nz_index] = variable_offset + i;
      ++nz_index;

      // theta0
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 0;
      ++nz_index;

      // kappa0
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 1;
      ++nz_index;

      // dkappa0
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 2;
      ++nz_index;

      // y0
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 4;
      ++nz_index;

      // theta1
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 5;
      ++nz_index;

      // kappa1
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 6;
      ++nz_index;

      // dkappa1
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 7;
      ++nz_index;

      // y1
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_index + 9;
      ++nz_index;

      // s
      iRow[nz_index] = i * 2 + 1;
      jCol[nz_index] = variable_offset + i;
      ++nz_index;
    }

    int constraint_offset = 2 * (SourceData.number - 1);
    for (int i = 0; i < SourceData.number; ++i) {
      iRow[nz_index] = constraint_offset + i;
      jCol[nz_index] = i * 5 + 3;
      ++nz_index;

      iRow[nz_index] = constraint_offset + i;
      jCol[nz_index] = i * 5 + 4;
      ++nz_index;
    }
  } else {
    std::fill(values, values + nele_jac, 0.0);
    // first, positional equality constraints
    int nz_index = 0;

    for (int i = 0; i + 1 < SourceData.number; ++i) {
      double delta_si = x[5 * SourceData.number + i];
      double theta_i = x[5 * i + 0];
      double dtheta_i = x[5 * i + 1];
      double ddtheta_i = x[5 * i + 2];
      double theta_i_1 = x[5 * (i + 1) + 0];
      double dtheta_i_1 = x[5 * (i + 1) + 1];
      double ddtheta_i_1 = x[5 * (i + 1) + 2];
      CoefInitialize(theta_i, dtheta_i, ddtheta_i, theta_i_1, dtheta_i_1,
                     ddtheta_i_1, delta_si);
      CoefDeriveInitialize(theta_i, dtheta_i, ddtheta_i, theta_i_1, dtheta_i_1,
                           ddtheta_i_1, delta_si);

      //这一部分是求gauss-legendre数值积分
      vector<double> epsiloni = {-0.90617984594, -0.53846931011, 0,
                                 0.53846931011, 0.90617984594};
      vector<double> wi = {0.236926885, 0.4786286705, 0.56888888889,
                           0.4786286705, 0.236926885};
      double gauss_legendre_x = 0;
      double gauss_legendre_y = 0;
      for (int j = 0; j < 5; ++j) {
        double segment_s = 0.5 * delta_si * epsiloni[j] + 0.5 * delta_si;
        double theta_segment_s = SpiralCurve(0, 1, segment_s);
        gauss_legendre_x += wi[j] * cos(theta_segment_s);
        gauss_legendre_y += wi[j] * sin(theta_segment_s);
      }
      gauss_legendre_x *= delta_si / 2;
      gauss_legendre_y *= delta_si / 2;

      auto posx_theta0 = GaussLegendreXDerivative(0, SourceData.m, delta_si);
      auto posx_kappa0 = GaussLegendreXDerivative(1, SourceData.m, delta_si);
      auto posx_dkappa0 = GaussLegendreXDerivative(2, SourceData.m, delta_si);
      auto posx_theta1 = GaussLegendreXDerivative(3, SourceData.m, delta_si);
      auto posx_kappa1 = GaussLegendreXDerivative(4, SourceData.m, delta_si);
      auto posx_dkappa1 = GaussLegendreXDerivative(5, SourceData.m, delta_si);
      auto posx_delta_s = GaussLegendreXDerivative(6, SourceData.m, delta_si);

      auto posy_theta0 = GaussLegendreYDerivative(0, SourceData.m, delta_si);
      auto posy_kappa0 = GaussLegendreYDerivative(1, SourceData.m, delta_si);
      auto posy_dkappa0 = GaussLegendreYDerivative(2, SourceData.m, delta_si);
      auto posy_theta1 = GaussLegendreYDerivative(3, SourceData.m, delta_si);
      auto posy_kappa1 = GaussLegendreYDerivative(4, SourceData.m, delta_si);
      auto posy_dkappa1 = GaussLegendreYDerivative(5, SourceData.m, delta_si);
      auto posy_delta_s = GaussLegendreYDerivative(6, SourceData.m, delta_si);

      // theta0
      values[nz_index] += (-posx_theta0);
      ++nz_index;

      // kappa0
      values[nz_index] += (-posx_kappa0);
      ++nz_index;

      // dkappa0
      values[nz_index] += (-posx_dkappa0);
      ++nz_index;

      // x0
      values[nz_index] += (-1.0);
      ++nz_index;

      // theta1
      values[nz_index] += (-posx_theta1);
      ++nz_index;

      // kappa1
      values[nz_index] += (-posx_kappa1);
      ++nz_index;

      // dkappa1
      values[nz_index] += (-posx_dkappa1);
      ++nz_index;

      // x1
      values[nz_index] += (1.0);
      ++nz_index;

      // delta_s
      values[nz_index] += (-posx_delta_s);
      ++nz_index;

      // for y coordinate
      // theta0
      values[nz_index] += (-posy_theta0);
      ++nz_index;

      // kappa0
      values[nz_index] += (-posy_kappa0);
      ++nz_index;

      // dkappa0
      values[nz_index] += (-posy_dkappa0);
      ++nz_index;

      // y0
      values[nz_index] += (-1.0);
      ++nz_index;

      // theta1
      values[nz_index] += (-posy_theta1);
      ++nz_index;

      // kappa1
      values[nz_index] += (-posy_kappa1);
      ++nz_index;

      // dkappa1
      values[nz_index] += (-posy_dkappa1);
      ++nz_index;

      // y1
      values[nz_index] += (1.0);
      ++nz_index;

      // delta_s
      values[nz_index] += (-posy_delta_s);
      ++nz_index;
    }

    for (int i = 0; i < SourceData.number; ++i) {
      values[nz_index] =
          2.0 * (x[i * 5 + 3] - SourceData.input_points[i].first);
      ++nz_index;

      values[nz_index] =
          2.0 * (x[i * 5 + 4] - SourceData.input_points[i].second);
      ++nz_index;
    }
  }
  return true;
}

// return the structure or values of the hessian
bool SpiralSmootherNLP::eval_h(Ipopt::Index n, const Ipopt::Number* x,
                               bool new_x, Ipopt::Number obj_factor,
                               Ipopt::Index m, const Ipopt::Number* lambda,
                               bool new_lambda, Ipopt::Index nele_hess,
                               Ipopt::Index* iRow, Ipopt::Index* jCol,
                               Ipopt::Number* values) {
  return true;
}

void SpiralSmootherNLP::finalize_solution(
    Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
    const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
    const Ipopt::Number* g, const Ipopt::Number* lambda,
    Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
    Ipopt::IpoptCalculatedQuantities* ip_cq) {
  // here is where we would store the solution to variables, or write to a file,
  // etc so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl
            << std::endl
            << "Solution of the primal variables, x" << std::endl;
  for (Ipopt::Index i = 0; i < SourceData.number; i++) {
    opt_thetai.emplace_back(x[5 * i]);
    opt_dthetai.emplace_back(x[5 * i + 1]);
    opt_ddthetai.emplace_back(x[5 * i + 2]);
    opt_xi.emplace_back(x[5 * i + 3]);
    opt_yi.emplace_back(x[5 * i + 4]);
  }
  for (Ipopt::Index i = 0; i < SourceData.number - 1; i++) {
    opt_deltas.emplace_back(x[5 * SourceData.number + i]);
  }
  vector<double> X;
  vector<double> Y;

  vector<double> X2;
  vector<double> Y2;

  for (auto& p : SourceData.input_points) {
    X2.emplace_back(p.first);
    Y2.emplace_back(p.second);
  }
  for (int i = 0; i < SourceData.number - 1; i++) {
    CoefInitialize(opt_thetai[i], opt_dthetai[i], opt_ddthetai[i],
                   opt_thetai[i + 1], opt_dthetai[i + 1], opt_ddthetai[i + 1],
                   opt_deltas[i]);
    for (int p = 0; p < 10; p++) {
      double delta_ssi = p * 0.1 * opt_deltas[i];

      //这一部分是求gauss-legendre数值积分
      vector<double> epsiloni = {-0.90617984594, -0.53846931011, 0,
                                 0.53846931011, 0.90617984594};
      vector<double> wi = {0.236926885, 0.4786286705, 0.56888888889,
                           0.4786286705, 0.236926885};
      double gauss_legendre_x = 0;
      double gauss_legendre_y = 0;
      for (int j = 0; j < SourceData.m; ++j) {
        double segment_s = 0.5 * delta_ssi * epsiloni[j] + 0.5 * delta_ssi;
        double theta_segment_s = SpiralCurve(0, 1, segment_s);
        gauss_legendre_x += wi[j] * cos(theta_segment_s);
        gauss_legendre_y += wi[j] * sin(theta_segment_s);
      }
      gauss_legendre_x *= delta_ssi / 2;
      gauss_legendre_y *= delta_ssi / 2;
      //

      X.push_back(opt_xi[i] + gauss_legendre_x);
      Y.push_back(opt_yi[i] + gauss_legendre_y);

      plt::clf();
      std::map<std::string, std::string> keywords;
      keywords.insert(pair<string, string>("color", "r"));
      keywords.insert(pair<string, string>("label", "Spiral Path"));
      plt::plot(X, Y, keywords);

      std::map<std::string, std::string> keywords1;
      keywords1.insert(pair<string, string>("label", "Pramal Points"));
      keywords1.insert(pair<string, string>("marker", "p"));
      plt::scatter(X2, Y2, 10.0, keywords1);

      std::map<std::string, std::string> keywords2;
      keywords2.insert(pair<string, string>("label", "Spiral Points"));
      keywords2.insert(pair<string, string>("marker", "^"));
      plt::scatter(opt_xi, opt_yi, 10.0, keywords2);

      plt::legend();
      plt::pause(0.05);
    }
  }

  // std::cout << std::endl
  //           << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  // for (Ipopt::Index i = 0; i < n; i++) {
  //   std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  // }
  // for (Ipopt::Index i = 0; i < n; i++) {
  //   std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  // }

  //打印目标变量的最终结果
  for (Index i = 0; i < n; i++) {
    std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  //打印各个约束条件的值
  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Ipopt::Index i = 0; i < m; i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }

  //打印最终的目标函数值
  cout << "f(x*) = " << obj_value << endl;
}
