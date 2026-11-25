/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.5.1 (Java 11 port): a tool for automated trimming of multiple
    sequence alignments, together with its readAl and statal companions.

    This file is part of the trimAl Java port, a derivative work of trimAl
    (https://github.com/inab/trimal).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

#include <iostream>
#include <string>
#include "alignment.h"
#include "sequencesMatrix.h"

int main(int argc, char **argv) {
  std::string mode(argv[1]);
  if (mode == "identity") {
    std::cout.precision(9);
    for (int i = 2; i < argc; i++) {
      alignment a;
      if (!a.loadAlignment(argv[i])) return 1;
      a.printColumnsIdentity_DescriptiveStats();
    }
  } else if (mode == "matrix") {
    std::string sequences[] = {"-A-C", "AC--", "--AC"};
    std::string names[] = {"alpha", "beta", "gamma"};
    sequencesMatrix matrix(sequences, names, 3, 4);
    int order[] = {2, 0, 1};
    int row[4];
    matrix.setOrder(order);
    for (int i = 0; i < 3; i++) {
      if (!matrix.getSequence(names[i], row)) return 1;
      std::cout << names[i];
      for (int j = 0; j < 4; j++) std::cout << " " << row[j];
      std::cout << std::endl;
    }
    int duplicateOrder[] = {1, 1, 2};
    matrix.setOrder(duplicateOrder);
    for (int i = 0; i < 3; i++) {
      matrix.getSequence(names[i], row);
      std::cout << names[i];
      for (int j = 0; j < 4; j++) std::cout << " " << row[j];
      std::cout << std::endl;
    }
  } else if (mode == "locale") {
    alignment a;
    std::cout << (a.loadAlignment(argv[2]) ? "loaded" : "not loaded") << std::endl;
    a.setOutputFormat(8, false);
    a.printAlignment();
  }
  return 0;
}
