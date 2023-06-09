/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-20 22:09:41
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-19 17:04:32
 * @FilePath: /SpiralSmoother/include/ipopt_interface.hpp
 * @Description: 螺旋线问题的非线性ipopt建模
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#ifndef IPOPT_INTERFACE_HPP
#define IPOPT_INTERFACE_HPP

#include <IpTNLP.hpp>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace Ipopt;
using namespace std;

// source waypoint data
class CsvSourceData {
 public:
  CsvSourceData(const char* filename) {
    GetSourceData(filename);
    cout << "等待平滑的路径点个数为：" << number << endl;
  }
  ~CsvSourceData() = default;

 private:
  void GetSourceData(const char* filename);

 public:
  int number;  //需要处理的点的数量
  vector<pair<double, double>> input_points;
  vector<double> distance;         //两点的距离
  vector<double> relative_theta_;  //存储原始点的heading
  double boundary = 0.4;           //位置约束
  int m = 5;                       //高斯-勒让德数值积分分段数量
};

// ipopt
class SpiralSmootherNLP : public Ipopt::TNLP {
 public:
  /** default constructor */
  SpiralSmootherNLP();

  /** default destructor */
  virtual ~SpiralSmootherNLP();

  CsvSourceData SourceData{"../source/path.csv"};

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m,
                            Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag,
                            Ipopt::TNLP::IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l,
                               Ipopt::Number* x_u, Ipopt::Index m,
                               Ipopt::Number* g_l, Ipopt::Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                  bool init_z, Ipopt::Number* z_L,
                                  Ipopt::Number* z_U, Ipopt::Index m,
                                  bool init_lambda, Ipopt::Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                           Ipopt::Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Index m, Ipopt::Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                          Ipopt::Index m, Ipopt::Index nele_jac,
                          Ipopt::Index* iRow, Ipopt::Index* jCol,
                          Ipopt::Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                      Ipopt::Number obj_factor, Ipopt::Index m,
                      const Ipopt::Number* lambda, bool new_lambda,
                      Ipopt::Index nele_hess, Ipopt::Index* iRow,
                      Ipopt::Index* jCol, Ipopt::Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can
   * store/write the solution */
  virtual void finalize_solution(
      Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
      const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
      const Ipopt::Number* g, const Ipopt::Number* lambda,
      Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq);
  //@}

 private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  PracticeNLP();
  SpiralSmootherNLP(const SpiralSmootherNLP&);
  SpiralSmootherNLP& operator=(const SpiralSmootherNLP&);
  //@}
};

#endif