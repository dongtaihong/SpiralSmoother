/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-20 22:18:16
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2023-04-19 16:47:33
 * @FilePath: /SpiralSmoother/src/SpiralSmoother.cpp
 * @Description: 螺旋线非线性问题的具体求解（问题建模由ipopt_interface完成）
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include "SpiralSmoother.hpp"

SpiralSmoother::SpiralSmoother() {
  Ipopt::SmartPtr<Ipopt::TNLP> spiral_smoother_nlp = new SpiralSmootherNLP();

  Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();

  app->Options()->SetStringValue("hessian_approximation", "limited-memory");
  app->Options()->SetNumericValue("tol", 1);
  app->Options()->SetNumericValue("acceptable_tol", 1);

  Ipopt::ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Ipopt::Solve_Succeeded) {
    std::cout << std::endl
              << std::endl
              << "*** Error during initialization!" << std::endl;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(spiral_smoother_nlp);

  if (status == Ipopt::Solve_Succeeded) {
    std::cout << std::endl
              << std::endl
              << "*** The problem solved!" << std::endl;
  } else {
    std::cout << std::endl
              << std::endl
              << "*** The problem FAILED!" << std::endl;
  }
}