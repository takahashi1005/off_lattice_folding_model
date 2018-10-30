#pragma once
#include <cmath>
#include <iostream>
#include <vector>

// 
template<typename T>
std::ostream& operator<< (std::ostream &os, std::vector<T> vec){
  for(int index=0;index<vec.size();++index){
    os << vec[index] << " ";
  }
  return os;
}

// 
template<typename T>
std::ostream& operator<< (std::ostream &os, std::vector<std::vector<T>> matrix){
  for(auto itr = matrix.begin();itr != matrix.end();++itr){
    for(auto jtr = itr->begin();jtr != itr->end();++jtr){
      os << *jtr << " ";
    }
    os << std::endl;
  }
  return os;
}

// ベクトル足し算
template<typename T>
std::vector<T> operator+ (std::vector<T> &left, std::vector<T> &right){
  std::vector<T> result(left.size());
  if(left.size() == right.size()){
    for(int index=0;index<left.size();++index){
      result[index] = left[index] + right[index];
    }
  } else {
   std::cout << "Wrong calculation, sum of vector\n";
   std::cout << left << std::endl;
   std::cout << right << std::endl;
  }
  return result;
}

//  ベクトル引き算
template<typename T>
std::vector<T> operator- (std::vector<T> &left, std::vector<T> &right){
  std::vector<T> result(left.size());
  if(left.size() == right.size()){
    for(int index=0;index<left.size();++index){
      result[index] = left[index] - right[index];
    }
  } else {
    std::cout << "Wrong calculation, minus of vector\n";
   std::cout << left << std::endl;
   std::cout << right << std::endl;
  }
  return result;
}

// ベクトルのスカラー倍
template<typename T>
std::vector<T> operator* (T coef, std::vector<T> &vec){
  std::vector<T> result(vec.size());
  for(int index=0;index<vec.size();++index){
    result[index] = coef * vec[index];
  }
  return result;
}
 
// ベクトルのスカラー倍（割り算）
template<typename T>
std::vector<T> operator/ (std::vector<T> &vec, T x){
  std::vector<T> result(vec.size());
  for(int index=0;index<vec.size();++index){
    result[index] = vec[index] / x;
  }
  return result;
}

// 内積
template<typename T>
T IP(std::vector<T> &left, std::vector<T> &right){
  T result = 0;
  if(left.size() == right.size()){
    for(int index=0;index<left.size();++index){
      result += left[index] * right[index];
    }
  } else {
    std::cout << "Wring calculation, inner product of vector\n";
   std::cout << left << std::endl;
   std::cout << right << std::endl;
  }
  return result;
}

// 外積
template<typename T>
std::vector<T> OP(std::vector<T> &left, std::vector<T> &right){
  std::vector<T> result(3,0);
    result[0] = left[1] * right[2] - left[2] * right[1];
    result[1] = left[2] * right[0] - left[0] * right[2];
    result[2] = left[0] * right[1] - left[1] * right[0];
  return result;
}

// 絶対値
template<typename T>
double abs(std::vector<T> vec){
  T result=0;
  for(int index=0;index<vec.size();++index){
    result += (double)(vec[index] * vec[index]);
  }
  return sqrt(result);
}

// なす角
template<typename T>
double angle(std::vector<T> &vec1, std::vector<T> &vec2){
  double theta;
  double cos;
  cos = IP(vec1,vec2)/abs(vec1)/abs(vec2);
  theta = acos(cos);
  return theta;
}
