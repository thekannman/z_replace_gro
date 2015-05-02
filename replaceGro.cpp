//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

#include "z_sim_params.hpp"
#include "z_string.hpp"
#include "z_vec.hpp"
#include "z_conversions.hpp"
#include "z_molecule.hpp"
#include "z_atom_group.hpp"
#include "z_gromacs.hpp"
#include "xdrfile_trr.h"
#include "boost/program_options.hpp"
#include <random>

using namespace arma;
// Units are nm, ps.

int main (int argc, char *argv[]) {
  po::options_description desc("Options");
  desc.add_options()
    ("help,h",  "Print help message and exit")
    ("group,g", po::value<std::string>()->default_value("SOL"),
     "Group from which molecules will be replaced")
    ("top", po::value<std::string>()->default_value("topol.top"),
     ".top file containing atomic/molecular properties")
    ("in,i", po::value<std::string>()->default_value("out.gro"),
     ".gro file containing list of atoms/molecules");
    ("out,o", po::value<std::string>()->default_value("conf.gro"),
     ".gro file to which output is written");
    ("description,d", po::value<std::string>()->default_value(""),
     "Description to include in output file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  std::vector<Molecule> molecules =
    gen_molecules(vm["top"].as<std::string>(), params);
  Atom_group all_atoms(vm["in"].as<std::string>(), molecules);
  Atom_group replaced_group;
  replaced_group.SetIndices(vm["group"].as<std::string>(),
                             select_group(groups,
                                          vm["group"].as<std::string>()),
                             all_atoms);
  std::string box_line = ReadLastLine(vm["gro"].as<std::string>());
  std::vector<std::string> box_line_split = split(box_line, " ");
  rowvec box;
  for (int i = 0; i < 3; i++)
    box << boost::lexical_cast<double>(split_line[i]);
  std::default_random_engine generator;
  std::uniform_int_distribution<unsigned> distribution(0,liquid_group.size()-1);
  double z_lower_bound = 1.00, z_upper_bound = 1.70;
  unsigned attempt_count=0;
  int rep_count = 0;
  const int kMaxReps = 10;
  const unsigned kMaxAttempts = 1e5;
  do {
    std::uniform_int_distribution<unsigned>
    distribution(0,replaced_group.size()-1);
    unsigned mol_to_replace i_mol = distribution(generator);
    if (replaced_group.positions(i_mol,2) > z_lower_bound &&
        replaced_group.positions(i_mol,2) < z_upper_bound &&
        replaced(i_mol) == 0) {
      //Insert repaced.group.replace_mol() call here.
      replaced(i) = 1;
      rep_count++;
    }
    attempt_count++;
  } while (rep_count < kMaxReps && attempt_count < kMaxAttempts);
  all_atoms.WriteGro(vm["out"].as<std::string>(), box,
                     vm["description"].as<std::string>());
} // main
