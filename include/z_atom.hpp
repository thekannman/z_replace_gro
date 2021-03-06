// Copyright (c) 2015 Zachary Kann
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// ---
// Author: Zachary Kann

// Contains information about a single atom. Typically used as part of a
// Molecule class.

#include <map>
#include <string>

#ifndef _Z_ATOM_HPP_
#define _Z_ATOM_HPP_

class Atom {
 public:
  // Better than the << operator for printing to std::cout.
  void Print() const;

  friend std::ostream &operator<<( std::ostream &output, const Atom &atom );

  // Basically a mutator for cross_sigma and cross_epsilon
  inline void set_cross_lj(const std::string& other_atom,
                           const double sigma, const double epsilon) {
    cross_sigma_[other_atom] = sigma;
    cross_epsilon_[other_atom] = epsilon;
  }

  //Mutators
  inline void set_type(const std::string& type) { type_ = type; }
  inline void set_name(const std::string& name) { name_ = name; }
  inline void set_mass(const double mass) { mass_ = mass; }
  inline void set_charge(const double charge) { charge_ = charge; }
  inline void set_sigma(const double sigma) { sigma_ = sigma; }
  inline void set_epsilon(const double epsilon) { epsilon_ = epsilon; }

  //Accessors
  inline std::string name() const { return name_; }
  inline std::string type() const { return type_; }
  inline double mass() const { return mass_; }
  inline double charge() const { return charge_; }
  inline double sigma() const { return sigma_; }
  inline double epsilon() const { return epsilon_; }

 private:
  std::string name_;
  std::string type_;
  double mass_;
  double charge_;
  double sigma_;
  double epsilon_;
  std::map<std::string, double> cross_sigma_;
  std::map<std::string, double> cross_epsilon_;
};

#endif
