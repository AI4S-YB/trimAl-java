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
