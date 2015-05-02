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
#include "z_file.hpp"
#include "z_molecule.hpp"
#include "z_atom_group.hpp"
#include "z_gromacs.hpp"
#include "boost/program_options.hpp"
#include <random>

namespace po = boost::program_options;
// Units are nm, ps.

int main (int argc, char *argv[]) {
  SimParams params;
  po::options_description desc("Options");
  desc.add_options()
    ("help,h",  "Print help message and exit")
    ("count,c", po::value<int>()->default_value(0),
     "Number of molecules to replace")
    ("group,g", po::value<std::string>()->default_value("SOL"),
     "Group from which molecules will be replaced")
    ("new,n", po::value<std::string>()->default_value("He"),
     "Name of new molecule.")
    ("top", po::value<std::string>()->default_value("topol.top"),
     ".top file containing atomic/molecular properties")
    ("index,n", po::value<std::string>()->default_value("index.ndx"),
     ".ndx file containing atomic indices for groups")
    ("in,i", po::value<std::string>()->default_value("out.gro"),
     ".gro file containing list of atoms/molecules")
    ("out,o", po::value<std::string>()->default_value("conf.gro"),
     ".gro file to which output is written")
    ("description,d", po::value<std::string>()->default_value(""),
     "Description to include in output file");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }
  std::map<std::string, std::vector<int> > groups;
  groups = ReadNdx(vm["index"].as<std::string>());

  std::vector<Molecule> molecules = GenMolecules(vm["top"].as<std::string>(),
                                                 params);
  AtomGroup all_atoms(vm["in"].as<std::string>(), molecules);
  Molecule* new_molecule;
  AtomGroup replaced_group(vm["group"].as<std::string>(),
                           SelectGroup(groups, vm["group"].as<std::string>()),
                           all_atoms);
  for (std::vector<Molecule>::iterator i_mol = molecules.begin();
       i_mol != molecules.end(); ++i_mol) {
    if ((*i_mol).name() == vm["new"].as<std::string>()) {
      new_molecule = &(*i_mol);
      break;
    }
    assert(i_mol+1 != molecules.end());
  }
  std::string box_line = ReadLastLine(vm["in"].as<std::string>());
  std::vector<std::string> box_line_split = Split(box_line, ' ');
  arma::rowvec box(DIMS);
  for (int i = 0; i < DIMS; i++)
    box(i) = boost::lexical_cast<double>(box_line_split[i]);
  std::default_random_engine generator;
  std::uniform_int_distribution<unsigned> distribution(0,
                                                       replaced_group.size()-1);
  //TODO(Zak): remove hard coding of bounds.
  double z_lower_bound = 1.00, z_upper_bound = 1.70;
  unsigned attempt_count = 0;
  int rep_count = 0;
  const int kMaxReps = vm["count"].as<int>();
  const unsigned kMaxAttempts = 1e5;
  all_atoms.UpdateCom();
  arma::irowvec replaced = arma::zeros<arma::irowvec>(replaced_group.size());
  do {
    std::uniform_int_distribution<unsigned>
    distribution(0,replaced_group.num_molecules()-1);
    unsigned i_mol = distribution(generator);
    if (replaced_group.position(i_mol,2) > z_lower_bound &&
        replaced_group.position(i_mol,2) < z_upper_bound &&
        replaced(i_mol) == 0) {
      all_atoms.ReplaceMolecule(*new_molecule, i_mol);
      replaced_group.RemoveMolecule(i_mol);
      replaced(i_mol) = 1;
      rep_count++;
    }
    attempt_count++;
  } while (rep_count < kMaxReps && attempt_count < kMaxAttempts);
  all_atoms.WriteGro(vm["out"].as<std::string>(), box,
                     vm["description"].as<std::string>());
} // main
