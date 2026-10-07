// Compile against the unmodified trimAl implementations; no alignment/main dependency.
#include <bits/stdc++.h>
#define private public
#include "statisticsConservation.h"
#undef private
#include "statisticsGaps.h"
#include "similarityMatrix.h"
#include "utils.h"

using namespace std;

static float asFloat(uint32_t x) { float v; memcpy(&v, &x, sizeof(v)); return v; }
static double asDouble(uint64_t x) { double v; memcpy(&v, &x, sizeof(v)); return v; }
static uint32_t bits(float v) { uint32_t x; memcpy(&x, &v, sizeof(x)); return x; }
static uint64_t bits(double v) { uint64_t x; memcpy(&x, &v, sizeof(x)); return x; }
static void floats(const string &label, const float *v, int n) {
  cout << '|' << label << '=';
  for (int i = 0; i < n; ++i) { if (i) cout << ','; cout << bits(v[i]); }
}
static void ints(const string &label, const int *v, int n) {
  cout << '|' << label << '=';
  for (int i = 0; i < n; ++i) { if (i) cout << ','; cout << v[i]; }
}
static void matrix(similarityMatrix &m, int mode) {
  if (mode == 1) m.defaultNTSimMatrix();
  else if (mode == 3) m.defaultAASimMatrix();
  else if (mode == 4) m.defaultNTDegeneratedSimMatrix();
  else m.alternativeSimilarityMatrices(1, mode - 10);
}

int main() {
  ios::sync_with_stdio(false);
  string line;
  while (getline(cin, line)) {
    istringstream in(line);
    string cmd, id; in >> cmd >> id;
    cout << id;
    if (cmd == "M") {
      int mode; string alphabet, path; in >> mode >> alphabet >> path;
      similarityMatrix m;
      bool ok = true;
      if (mode == 0) { vector<char> fn(path.begin(), path.end()); fn.push_back(0); ok = m.loadSimMatrix(fn.data()); }
      else matrix(m, mode);
      cout << "|ok=" << ok;
      if (ok) {
        vector<float> d;
        for (char a : alphabet) for (char b : alphabet) d.push_back(m.getDistance(a, b));
        floats("dist", d.data(), d.size());
      }
    } else if (cmd == "A") {
      similarityMatrix m; m.alternativeSimilarityMatrices(1, AAType);
      struct Layout { int *hash; float **sim; float **dist; int size; } layout;
      static_assert(sizeof(layout) == sizeof(m), "Unexpected reference matrix object layout");
      memcpy(&layout, &m, sizeof(layout));
      cout << "|size=" << layout.size;
      vector<float> sim, dist;
      for (int i = 0; i < layout.size; ++i) for (int j = 0; j < layout.size; ++j) {
        sim.push_back(layout.sim[i][j]); dist.push_back(layout.dist[i][j]);
      }
      floats("sim", sim.data(), sim.size()); floats("dist", dist.data(), dist.size());
    } else if (cmd == "U") {
      int cols; uint32_t valueBits, thresholdBits; in >> cols >> valueBits >> thresholdBits;
      vector<string> rows(2, string(cols, 'A'));
      statisticsConservation c(rows.data(), 2, cols, DNAType);
      for (int i = 0; i < cols; ++i) c.getMdkwVector()[i] = asFloat(valueBits);
      cout << "|cut=" << bits(c.calcCutPoint(0, asFloat(thresholdBits)));
    } else if (cmd == "G") {
      int n; in >> n;
      vector<string> rows(n, "");
      statisticsGaps g(rows.data(), n, 0, DNAType);
      cout << "|slope=" << g.calcCutPoint2ndSlope();
    } else if (cmd == "S") {
      int datatype, n, cols, hw, mode, secondWindow;
      uint32_t minBits, simBits, gapBits;
      in >> datatype >> n >> cols >> hw >> mode >> minBits >> simBits >> gapBits >> secondWindow;
      vector<string> rows(n); for (auto &r : rows) in >> r;
      statisticsGaps g(rows.data(), n, cols, datatype);
      bool gw = g.applyWindow(hw);
      ints("gaps", g.getGapsWindow(), cols);
      cout << "|gw=" << gw << "|slope=" << g.calcCutPoint2ndSlope()
           << "|gapcut=" << bits(g.calcCutPoint(asFloat(minBits), asFloat(gapBits)));
      similarityMatrix m; matrix(m, mode);
      statisticsConservation c(rows.data(), n, cols, datatype);
      c.setSimilarityMatrix(&m);
      cout << "|vectors=" << c.calculateVectors(rows.data(), g.getGapsWindow());
      floats("q", c.Q, cols); floats("mdk", c.MDK, cols);
      vector<float> identity;
      for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) identity.push_back(c.matrixIdentity[i][j]);
      floats("identity", identity.data(), identity.size());
      cout << "|cw=" << c.applyWindow(hw);
      floats("window", c.getMdkwVector(), cols);
      cout << "|conscut=" << bits(c.calcCutPoint(asFloat(minBits), asFloat(simBits)));
      if (secondWindow >= 0) {
        cout << "|gw2=" << g.applyWindow(secondWindow);
        ints("gaps2", g.getGapsWindow(), cols);
        cout << "|cw2=" << c.applyWindow(secondWindow);
        floats("window2", c.getMdkwVector(), cols);
      }
    } else if (cmd == "F") {
      char type; uint64_t raw; int precision, fixed;
      in >> type >> raw >> precision >> fixed;
      ostringstream out; out.precision(precision);
      if (fixed) out << std::fixed;
      if (type == 'f') out << asFloat(raw); else out << asDouble(raw);
      cout << "|text=" << out.str();
    } else if (cmd == "T") {
      uint32_t raw; in >> raw;
      float v = asFloat(raw);
      cout << "|exp=" << bits(exp(v)) << "|log10=" << bits(log10(v))
           << "|sqrt=" << bits(sqrt(v))
           << "|expDouble=" << bits((float) exp((double) v))
           << "|log10Double=" << bits((float) log10((double) v))
           << "|floatOverloads=" << is_same<decltype(exp(v)), float>::value;
    } else if (cmd == "I") {
      string hex, s; in >> hex;
      for (size_t i = 0; i < hex.size(); i += 2) s.push_back(stoul(hex.substr(i, 2), nullptr, 16));
      cout << "|atoi=" << atoi(s.c_str()) << "|atof=" << bits(atof(s.c_str()))
           << "|number=" << utils::isNumber(const_cast<char *>(s.c_str()));
    } else if (cmd == "R") {
      uint64_t raw; in >> raw; double v = asDouble(raw);
      cout << "|inf=" << utils::roundToInf(v) << "|nearest=" << utils::roundInt(v) << "|sup=" << utils::roundToSup(v);
    } else if (cmd == "Q") {
      char type; int n; in >> type >> n;
      if (type == 'f') {
        vector<float> v(n); uint32_t raw; for (auto &x : v) { in >> raw; x = asFloat(raw); }
        utils::quicksort(v.data(), 0, n - 1); floats("sorted", v.data(), n);
      } else if (type == 'p') {
        vector<array<int, 2>> v(n); vector<int *> p(n);
        for (int i = 0; i < n; ++i) { in >> v[i][0]; v[i][1] = i; p[i] = v[i].data(); }
        utils::quicksort(p.data(), 0, n - 1);
        vector<int> indices; for (auto x : p) indices.push_back(x[1]); ints("sorted", indices.data(), n);
      } else {
        vector<int> v(n); for (auto &x : v) in >> x;
        utils::quicksort(v.data(), 0, n - 1); ints("sorted", v.data(), n);
      }
    } else {
      cerr << "Unknown command " << cmd << '\n'; return 2;
    }
    cout << '\n';
  }
}
