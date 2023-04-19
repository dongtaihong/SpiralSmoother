/*
 * @Author: 董泰宏 2396203400@qq.com
 * @Date: 2022-12-20 22:33:52
 * @LastEditors: 董泰宏 2396203400@qq.com
 * @LastEditTime: 2022-12-22 14:46:42
 * @FilePath: /SpiralSmoother/include/SpiralSmoother.hpp
 * @Description: 螺旋线非线性问题的具体求解（问题建模由ipopt_interface完成）
 * Copyright (c) 2022 by 董泰宏 email: 2396203400@qq.com, All Rights Reserved.
 */
#include <IpIpoptApplication.hpp>
#include <iostream>

#include "ipopt_interface.hpp"

class SpiralSmoother
{
public:
    SpiralSmoother();
    ~SpiralSmoother() = default;
    void SpiralToCartesian();
};