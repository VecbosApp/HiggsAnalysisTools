#ifndef KFACTOREVALUATOR_H
#define KFACTOREVALUATOR_H

class kFactorEvaluator {
public:
  kFactorEvaluator() {}
  virtual ~kFactorEvaluator() {}
  virtual float evaluate(float pT) {}
};

// Higgs -> WW at different masses
// k-factor as a function of pt(WW) for gg->h (mh = 130 GeV)
// by G.Davatz updated by droll (30/09/05)
class ggHiggs120Evaluator : public kFactorEvaluator {
public:
  ggHiggs120Evaluator() {}
  ~ggHiggs120Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs130Evaluator : public kFactorEvaluator {
public:
  ggHiggs130Evaluator() {}
  ~ggHiggs130Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs140Evaluator : public kFactorEvaluator {
public:
  ggHiggs140Evaluator() {}
  ~ggHiggs140Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs150Evaluator : public kFactorEvaluator {
public:
  ggHiggs150Evaluator() {}
  ~ggHiggs150Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs160Evaluator : public kFactorEvaluator {
public:
  ggHiggs160Evaluator() {}
  ~ggHiggs160Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs165Evaluator : public kFactorEvaluator {
public:
  ggHiggs165Evaluator() {}
  ~ggHiggs165Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs170Evaluator : public kFactorEvaluator {
public:
  ggHiggs170Evaluator() {}
  ~ggHiggs170Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs180Evaluator : public kFactorEvaluator {
public:
  ggHiggs180Evaluator() {}
  ~ggHiggs180Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs190Evaluator : public kFactorEvaluator {
public:
  ggHiggs190Evaluator() {}
  ~ggHiggs190Evaluator() {}
  float evaluate(float pT);
};

class ggHiggs200Evaluator : public kFactorEvaluator {
public:
  ggHiggs200Evaluator() {}
  ~ggHiggs200Evaluator() {}
  float evaluate(float pT);
};

// VBF kFactor does not depend by pT
class VBFHiggsEvaluator : public kFactorEvaluator {
public:
  VBFHiggsEvaluator() {}
  ~VBFHiggsEvaluator() {}
  float evaluate(float pT) {return 1.0;}
};

// default constant kFactor does not depend by pT
class flatHiggsEvaluator : public kFactorEvaluator {
public:
  flatHiggsEvaluator() {}
  ~flatHiggsEvaluator() {}
  float evaluate(float pT) {return 1.0;}
};


// WW non resonant 
// Fabian

class WWEvaluator : public kFactorEvaluator {
public:
  WWEvaluator() {}
  ~WWEvaluator() {}
  float evaluate(float pT);
};

#endif //kFactorEvaluator
