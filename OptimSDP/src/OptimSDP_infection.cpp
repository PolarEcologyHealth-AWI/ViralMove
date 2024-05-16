// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <vector>
#include <cassert> 
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// Namespace for Parameters
namespace sdp {

template <class T>

class checkedVector : public std::vector<T>
{
public:
  
  checkedVector()
  { }
  checkedVector(size_t n, const T& value = T())
    : std::vector<T>(n,value)
    { }
  checkedVector(T* i, T* j)
    : std::vector<T>(i,j)
    { }
  
  T& operator [] (ptrdiff_t index)
  {
    assert (index >= 0 && index < static_cast <ptrdiff_t> (size()));
    return std::vector<T>::operator[](index);
  }
};

// ---
// 2d
// ---

template <class T>
class Matrix
  : public sdp::checkedVector< sdp::checkedVector<T> >
{
protected:
  size_t rows,
  columns;
public:
  Matrix(size_t x = 0, size_t y = 0)
    : sdp::checkedVector< sdp::checkedVector<T> > (x,
      sdp::checkedVector<T>(y)), rows(x), columns(y)
      {}
  
  size_t Rows() const {return rows;}
  size_t Columns() const {return columns;}
  
  void init(const T& Value) {
    for (size_t i=0; i< rows; ++i)
      for (size_t j=0; j < columns; ++j)
        sdp::Matrix<T>::operator[](i)(j) = Value;
  }
  
  void resize(size_t x, size_t y, T t=T()) {
    sdp::checkedVector< sdp::checkedVector<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix<T>::operator[](i).resize(y, t);
    rows = x ; columns = y;
  }
  
  Matrix<T>& I()
  {
    for(size_t i = 0; i < rows; ++i)
      for(size_t j = 0; j < columns; ++j)
        sdp::Matrix<T>::operator[](i)(j) = (i == j) ? T(1) : T(0);
    return *this;
  }
  
};

// ---
// 3d
// ---

template <class T>
class Matrix3D : public sdp::checkedVector<Matrix<T> >
{
protected:
  size_t rows,
  columns,
  zDim;
public:
  Matrix3D(size_t x = 0, size_t y = 0, size_t z = 0)
    : sdp::checkedVector<Matrix<T> >(x, Matrix<T>(y, z)),
      rows(x), columns(y), zDim(z)
      {}
  
  size_t Rows()   const {return rows;}
  size_t Columns() const {return columns;}
  size_t zDIM()   const {return zDim;}
  
  void init(const T& value)
  {
    for (size_t i = 0; i < rows; ++i)
      sdp::Matrix3D<T>::operator[](i).init(value);
  }
  
  void resize (size_t x, size_t y, size_t z, T t=T()) {
    sdp::checkedVector< Matrix<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix3D<T>::operator[](i).resize(y,z,t);
    rows = x ; columns = y; zDim = z;
  }
};


// ---
// 4d
// ---

template <class T>
class Matrix4D : public sdp::checkedVector<Matrix3D<T> >
{
protected:
  size_t rows,
  columns,
  zDim,
  z2Dim;
public:
  Matrix4D(size_t x = 0, size_t y = 0, size_t z = 0, size_t z2 = 0)
    : sdp:: checkedVector<Matrix3D<T> >(x, Matrix3D<T>(y, z, z2)),
      rows(x), columns(y), zDim(z), z2Dim(z2)
      {}
  
  size_t Rows()    const {return rows;}
  size_t Columns() const {return columns;}
  size_t zDIM()    const {return zDim;}
  size_t z2DIM()   const {return z2Dim;}
  
  void init(const T& value)
  {
    for (size_t i = 0; i < rows; ++i)
      sdp::Matrix4D<T>::operator[](i).init(value);
  }
  
  void resize (size_t x, size_t y, size_t z, size_t z2, T t=T())
  {
    sdp::checkedVector< Matrix3D<T> >::resize(x);
    for (size_t i = 0; i < x; ++i)
      sdp::Matrix4D<T>::operator[](i).resize(y,z,z2,t);
    rows = x ; columns = y; zDim = z; z2Dim = z2;
  }
};




// Init
int direction;
int MinT;
int MaxT;
int NSites;
int MaxX;

// Individual
double B0;
double w;
double xc;
arma::vec expFor;
arma::vec expFly;
arma::vec b0;
arma::vec b1;
arma::vec b2;
double pred_a1;
double pred_a2;
double c;
double speed;
arma::vec WindAssist;
arma::vec WindProb;
arma::vec ZStdNorm;
arma::vec PStdNorm;
double decError;

// Sites
arma::mat dist;
arma::mat bear;
arma::vec nTR_x;
arma::vec nTR_y;
arma::mat expend;
arma::vec y_gain;
arma::vec penalty;

// Output
sdp::Matrix4D<float> FMatrix;
sdp::Matrix4D<float> DMatrix1;
sdp::Matrix4D<float> DMatrix2;
sdp::Matrix4D<float> PMatrix1;
sdp::Matrix4D<float> PMatrix2;

// Internal paramters
double Fintensity = 0.0;

} // end spd Namespace



// Eval (not exported):
// Function use to truncate value to range between min and max
int Eval(
    double x,
    double min,
    double max
)
{
  return (x <= min) ? 0 : ((x >= (max)) ? 2 : 1);
}


// Decision Error (not exported)
// Function to calculate decision error based on the importance of the decision
double DecError(
    double var,
    double decError
)
{
  double help1 = 1.0 + (decError * var);
  double help  = 1.0 / help1;
  return help;
}


// Round (not exported)
// Function equivalent to round(value, 0) in R
int Round(double value)
{
  if (value - floor(value) < 0.5) return (int)(value);
  else return ((int)(value) + 1);
}


// Sigmoidal TR function (not exported)
inline double FNtanh(double g)
{
  return ((( (exp(g) - exp(-g)) / (exp(g) + exp(-g)) ) + 1.0)/ 2.0);
}


// InterpolateTR
// Function to derive Terminal Reward for t
double InterpolateTR (int c_tr,
                      arma::vec  c_nTR_x,
                      arma::vec  c_nTR_y)
{
  int i = 0;
  double help = 0.0;
  double help_t = (double)(c_tr);
  double min = 0.0; double max = 2.0;
  
  while((c_nTR_x(i) < help_t) && (i < (c_nTR_x.size()))) ++i;
  if((i==0) || (c_nTR_x(i) - c_nTR_x[i-1]) == 0)
    help = c_nTR_y(i);
  else
    help = c_nTR_y(i-1) +
      (c_nTR_y(i) - c_nTR_y(i-1))/
        (c_nTR_x(i) - c_nTR_x(i-1)) * (help_t - c_nTR_x(i-1));
  
  if (help < min) help = min;
  if (help > max) help = max;
  return help;
}


// CalculateTR
// Function to calculate the Terminal Reward
double CalculateTR(int c_x,
                   int c_t,
                   double c_w,
                   double c_xc,
                   double c_B0,
                   arma::vec  c_nTR_x,
                   arma::vec  c_nTR_y)
{
  double TimeReward, StateReward, TR_tmp = 0.0;
  
  TimeReward = InterpolateTR(c_t, c_nTR_x, c_nTR_y);
  if (c_x == 0) StateReward = 0.0;
  else StateReward = FNtanh(c_w * ((double)(c_x) - c_xc));
  
  TR_tmp = TimeReward * StateReward + c_B0;
  return TR_tmp;
}

////////////////////////////////////////////////////////////////
////// Functions used in Backward and Forward Simulation ///////
////////////////////////////////////////////////////////////////


// Predation
// Calculates the reduction in fitness based on predation
double Predation (
    const int& time,
    const int& site,
    const int& x,
    const int& inf,
    double u,
    double f
)
{
  const double rew_tol = 0.0000005;
  double netgain, cor_u, cor_x, red;
  
  if (inf > 0) {
    red = 1 + sdp::expFor[0]/100;
  } else {
    red = 1;
  }
  
  cor_u = u;
  cor_x = (double)(x);
  netgain = ((u * f) - (sdp::expend(site, time)*red));
  if (netgain <= 0.0) netgain = rew_tol;
  if (cor_u   <= 0.0) cor_u   = rew_tol;
  if (cor_x   <= 0.0) cor_x   = rew_tol;
  
  double help1 = exp((sdp::pred_a1 + 1) * log(cor_x + netgain));
  double help2 = exp((sdp::pred_a1 + 1) * log(cor_x));
  double help3 = (sdp::pred_a1 + 1.0) * netgain;
  double help4 = exp(sdp::pred_a2 * log(cor_u));
  double help5 = sdp::b0(site) + sdp::b1(site) * ((help1 - help2)/help3) * sdp::b2(site) * help4;
  
  return help5;
} // end Predation


// FindFitness value
// Interpolates FitnessValue for time, site, x, inf and feeding intensity
double FindF (
    const int& time, 
    const int& site, 
    const int& x,
    const int& inf,
    int accuracy, 
    double gain, 
    double fx, 
    double u
)
{
  double res1, res2, part1, part2, interpolReward, red;
  res1 = 0.0; res2 = 0.0; part1 = 0.0; part2 = 0.0;
  
  if (inf > 0) {
    red = 1 + sdp::expFor[0]/100;
  } else {
    red = 1;
  }
  
  
  double expenditure = sdp::expend(site, time)*red;
  double nextx = (double)(x) + gain * u - expenditure;
  switch (Eval(nextx, 0.0, (double)(sdp::MaxX)))
  {
  case 0: interpolReward = sdp::FMatrix[time+1][site][0][inf]; break;
  case 1: {    res1 = nextx - (int)(nextx);
    res2 = (int)(nextx) + 1.0 - nextx;
    
    part1 = res1 * sdp::FMatrix[time+1][site][(int)(nextx)+1][inf];
    part2 = res2 * sdp::FMatrix[time+1][site][(int)(nextx)][inf];
    interpolReward = part1 + part2;
    break;
  }    
  case 2: interpolReward = sdp::FMatrix[time+1][site][sdp::MaxX][inf]; break;
  }
  double fx_new = fx + (sdp::PStdNorm(accuracy) * interpolReward);
  return fx_new;
} // end FindF


// Foraging
// Major Foraging function: Optimization of Foraging intensity
double Foraging(const int& time,
                const int& site,
                const int& x,
                const int& inf
)
{
  const double r = 0.61803399;
  const double tol = 0.0000005;
  double mean, gain, SD, Freward;
  double c, u0, u1, u2, u3, f1, f2, hold1, hold2;
  double hold1_old, hold2_old, f1_old, f2_old;
  int NStdNorm = sdp::ZStdNorm.size();
  
  c  = 1.0 - r;
  u0 = 0.0;
  u1 = r;
  u3 = 1.0;
  u2 = u1 + c * (u3 - u1);
  
  f1 = 0.0; f2 = 0.0;
  
  mean = sdp::y_gain(site);
  SD   = 1;
  
  for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
  {
    gain = mean + sdp::ZStdNorm(accuracy) * SD;
    f1_old = f1;
    f2_old = f2;
    f1 = FindF(time, site, x, inf, accuracy, gain, f1_old, u1);
    f2 = FindF(time, site, x, inf, accuracy, gain, f2_old, u2);
  }
  
  f1_old = f1;
  f2_old = f2;
  f1 = (1.0 - Predation(time, site, x, inf, u1, f1)) * f1;
  f2 = (1.0 - Predation(time, site, x, inf, u2, f2)) * f2;
  
  while ((fabs(u3 - u0)) > tol)
  {
    if (f2 > f1)
    {
      u0 = u1;
      u1 = u2;
      u2 = (r * u1) + (c * u3);
      hold2 = 0.0;
      for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
      {
        gain = mean+ sdp::ZStdNorm(accuracy) * SD;
        hold2_old = hold2;
        hold2 = FindF(time, site, x, inf, accuracy, gain, hold2_old, u2);
      }
      f1 = f2;
      f2 = ((1.0 - Predation(time, site, x, inf, u2, hold2)) * hold2);
    }
    if (f2 <= f1)
    {
      u3 = u2;
      u2 = u1;
      u1 = (r * u2) + (c * u0);
      hold1 = 0.0;
      for (int accuracy = 0; accuracy < NStdNorm; ++accuracy)
      {
        gain = mean + sdp::ZStdNorm(accuracy) * SD;
        hold1_old = hold1;
        hold1 = FindF(time, site, x, inf, accuracy, gain, hold1_old, u1);
      }
      f2 = f1;
      f1 = ((1.0 - Predation(time, site, x, inf, u1, hold1)) * hold1);
    }
  }  //end of while loop
  
  if (f1 <= f2)
  {
    sdp::Fintensity = u2;
    Freward = f2 * sdp::penalty(site);
  }
  else //(f1 > f2)
  {
    sdp::Fintensity = u1;
    Freward = f1 * sdp::penalty(site);
  }
  return Freward;
} // end Foraging

// Flying
double Flying(
    const int& time,
    const int& site,
    const int& x,
    const int& inf,
    int dep_site)
{
  double totalD, Whold, range, Sqr_c, Sqr_ca;
  double distance, t, interpolReward, nextx, red;
  double res1, res2, tim1, tim2, part1, part2, part3, part4;
  
  if (inf > 0) {
    red = 1 - sdp::expFly[0]/100;
  } else {
    red = 1;
  }
  
  
  
  interpolReward = 0.0;
  totalD = 0.0;
  Whold = 0.0;
  
  range = (sdp::c * red) * (1.0 - (1.0/ (sqrt(1.0 + ( ((double)(x)) / (double)(sdp::MaxX))))));
  
  totalD = sdp::dist(site, dep_site);
  for (int h = 0; h < sdp::WindProb.size(); ++h)
  {
    distance = totalD * (1.0 + sdp::WindAssist(h));
    
    Sqr_c  = ((sdp::c * red) * (sdp::c *red));
    Sqr_ca = ((sdp::c * red) - (range - distance))*((sdp::c * red) - (range - distance));
    nextx = ((Sqr_c/Sqr_ca) - 1.0) * sdp::MaxX;
    t = (double)(time) + (distance / sdp::speed);
    
    if (t >= (sdp::MaxT-sdp::MinT)) interpolReward = 0.0;
    else  if (nextx <= 0.0) interpolReward = 0.0;
    else
    {
      res1 = nextx - (int)(nextx);
      res2 = (int)(nextx) + 1.0 - nextx;
      tim1 = t - (int)(t);
      tim2 = (int)(t) + 1.0 - t;
      
      if ((nextx+1.0) > (double)(sdp::MaxX))
      {
        part1 = tim1*res1*sdp::FMatrix[(int)(t)+1][dep_site][sdp::MaxX][inf];
        part3 = tim2*res1*sdp::FMatrix[(int)(t)][dep_site][sdp::MaxX][inf];
      }
      else
      {
        part1 = tim1*res1*sdp::FMatrix[(int)(t)+1][dep_site][(int)(nextx)+1][inf];
        part3 = tim2*res1*sdp::FMatrix[(int)(t)][dep_site][(int)(nextx)+1][inf];
      }
      part2 = tim1*res2*sdp::FMatrix[(int)(t)+1][dep_site][(int)(nextx)][inf];
      part4 = tim2*res2*sdp::FMatrix[(int)(t)][dep_site][(int)(nextx)][inf];
      interpolReward = part1+part2+part3+part4;
    }
    double help_old = Whold;
    Whold = help_old + sdp::WindProb(h) * interpolReward;
  } //End h
  return Whold;
}  // end Flying



////////////////////////////////////////////////////////////////
////// Export Functions: Backward Iteration ////////////////////
////////////////////////////////////////////////////////////////


// [[Rcpp::export]]
void Init( int direction,
           int MinT,
           int MaxT,
           int NSites,
           int MaxX,
           double w,
           double xc,
           double expFor,
           double expFly,
           double B0,
           Rcpp::NumericVector b0,
           Rcpp::NumericVector b1,
           Rcpp::NumericVector b2,
           double pred_a1,
           double pred_a2,
           double c,
           double speed,
           Rcpp::NumericVector WindAssist,
           Rcpp::NumericVector WindProb,
           Rcpp::NumericVector ZStdNorm,
           Rcpp::NumericVector PStdNorm,
           Rcpp::NumericVector nTR_x,
           Rcpp::NumericVector nTR_y,
           double decError,
           arma::mat dist,
           arma::mat bear,
           Rcpp::NumericVector y_gain,
           arma::mat y_expend,
           Rcpp::NumericVector penalty
)
{
  // Basic Parms
  sdp::direction = direction;
  sdp::MinT   = MinT;
  sdp::MaxT   = MaxT;
  sdp::NSites = NSites;
  sdp::MaxX   = MaxX;
  
  // Individual
  sdp::B0 = B0;
  sdp::w = w;
  sdp::xc = xc;
  sdp::expFor = expFor;
  sdp::expFly = expFly;
  sdp::b0 = b0;
  sdp::b1 = b1;
  sdp::b2 = b2;
  sdp::pred_a1 = pred_a1;
  sdp::pred_a2 = pred_a2;
  sdp::c = c;
  sdp::speed = speed;
  sdp::WindAssist = WindAssist;
  sdp::WindProb = WindProb;
  sdp::ZStdNorm = ZStdNorm;
  sdp::PStdNorm = PStdNorm;
  sdp::nTR_x = nTR_x;
  sdp::nTR_y = nTR_y;
  sdp::decError = decError;
  
  // Sites
  sdp::dist   = dist;
  sdp::bear   = bear;
  sdp::y_gain = y_gain;
  sdp::expend = y_expend;
  sdp::penalty = penalty;
  
  sdp::Fintensity = 0.0;
}


// [[Rcpp::export]]
Rcpp::List BackwardIteration() {
  
  /// Terminal Reward
  sdp::Matrix4D<float>FM_TR;
  FM_TR.resize((sdp::MaxT-sdp::MinT)+1, sdp::NSites+1, sdp::MaxX+1, 2, 0.0);
  
  
  for (int r = 0; r <= (sdp::MaxT-sdp::MinT); ++r)
  {
    for (int c = 0; c <= sdp::MaxX; ++c)
    {
      FM_TR[r][sdp::NSites][c][0] = CalculateTR(c, r, sdp::w, sdp::xc, sdp::B0, sdp::nTR_x, sdp::nTR_y);
      FM_TR[r][sdp::NSites][c][1] = CalculateTR(c, r, sdp::w, sdp::xc, sdp::B0, sdp::nTR_x, sdp::nTR_y);
    }
  }

  for (int r = 0; r < sdp::NSites; ++r)
  {
    for (int c = 0; c <= sdp::MaxX; ++c)
    {
      FM_TR[(sdp::MaxT-sdp::MinT)][r][c][0] = sdp::B0;
      FM_TR[(sdp::MaxT-sdp::MinT)][r][c][1] = sdp::B0;
    }
  }
  
  sdp::FMatrix  = FM_TR;
  
  sdp::Matrix4D<float>DMatrix1;
  sdp::DMatrix1.resize(sdp::NSites+1, (sdp::MaxT-sdp::MinT)+1,sdp::MaxX+1, 2, 0.0);
  sdp::Matrix4D<float>DMatrix2;
  sdp::DMatrix2.resize(sdp::NSites+1, (sdp::MaxT-sdp::MinT)+1,sdp::MaxX+1, 2, 0.0);
  sdp::Matrix4D<float>PMatrix1;
  sdp::PMatrix1.resize(sdp::NSites+1, (sdp::MaxT-sdp::MinT)+1,sdp::MaxX+1, 2, 0.0);
  sdp::Matrix4D<float>PMatrix2;
  sdp::PMatrix2.resize(sdp::NSites+1, (sdp::MaxT-sdp::MinT)+1,sdp::MaxX+1, 2, 0.0);
  
  
  Rcpp::List out(5);
  
  arma::vec  Mrew(sdp::NSites+1);
  double     max_Reward, decision;

  
  // // Progress bar
  // Progress p(sdp::MaxT*(sdp::NSites), pbar);
  
  for (int inf = 0; inf <=1; ++inf)
  {
    for (int time = ((sdp::MaxT-sdp::MinT)-1); time >= 0; --time)
    {
      for (int site = 0; site < (sdp::NSites); ++site)
      {
  
        sdp::FMatrix[time][site][0][inf] = 0.0;
        // p.increment(); // update progress
  
        for(int x = 1; x <= (sdp::MaxX); ++x) {
  
          sdp::Fintensity = 0.0;
          max_Reward = 0.0;
          decision = 0.0;
          
          
          double f_help = Foraging(time, site, x, inf);
          sdp::FMatrix[time][site][x][inf] = f_help;
  
          decision   = - sdp::Fintensity - 1.0;
          sdp::DMatrix1[site][time][x][inf] = decision;
  
          max_Reward = f_help;
  
          for (int bb = 0; bb <= sdp::NSites; ++bb) Mrew(bb) = 0.0;
          
            for (int dest = site; dest <= sdp::NSites; ++dest)
            {
              ///if bearing TO BE INCLUDED
              float help_z = Flying(time, site, x, inf, dest);
              if(dest == site) help_z = 0;
              Mrew[dest] = help_z;
            }
            
            double help_flying = 0.0;
            int best_bb = 0;
            for (int bb = 0; bb <= sdp::NSites; ++bb)
            {
              if (Mrew[bb] > help_flying)
              {
                help_flying = Mrew(bb);
                best_bb = bb;
                sdp::DMatrix2[site][time][x][inf] = (double)(bb);
                if (help_flying > max_Reward) max_Reward = help_flying;
              }
            }
            
            sdp::PMatrix1[site][time][x][inf] = DecError((max_Reward - sdp::FMatrix[time][site][x][inf])/max_Reward, sdp::decError);
            sdp::PMatrix2[site][time][x][inf] = DecError((max_Reward - Mrew(best_bb))/max_Reward, sdp::decError);
            
            double sum_action = sdp::PMatrix1[site][time][x][inf] + sdp::PMatrix2[site][time][x][inf];
            //
            sdp::PMatrix1[site][time][x][inf] /= sum_action;
            sdp::PMatrix2[site][time][x][inf] /= sum_action;
            
            sdp::FMatrix[time][site][x][inf]  = max_Reward;
     
        } // end x
       
      } // end site
    } // end time
  }
  
  out[0] = sdp::FMatrix;
  out[1] = sdp::DMatrix1;
  out[2] = sdp::DMatrix2;
  out[3] = sdp::PMatrix1;
  out[4] = sdp::PMatrix2;
  
  return out;
}



// ///////////////////////////////////////////////////////////////
// ////// Export Functions: Forward Simulation ///////////////////
// ///////////////////////////////////////////////////////////////



// [[Rcpp::export]]
void InitSim (int MinT,
              int MaxT,
              int NSites,
              int MaxX,
              double w,
              double xc,
              double expFly,
              double expFor,
              double B0,
              Rcpp::NumericVector b0,
              Rcpp::NumericVector b1,
              Rcpp::NumericVector b2,
              double pred_a1,
              double pred_a2,
              double c,
              double speed,
              Rcpp::NumericVector WindAssist,
              Rcpp::NumericVector WindProb,
              Rcpp::NumericVector ZStdNorm,
              Rcpp::NumericVector PStdNorm,
              Rcpp::NumericVector nTR_x,
              Rcpp::NumericVector nTR_y,
              double decError,
              arma::mat dist,
              arma::mat bear,
              Rcpp::NumericVector y_gain,
              arma::mat y_expend)
{
  // Basic Parms
  sdp::MinT   = MinT;
  sdp::MaxT   = MaxT;
  sdp::NSites = NSites;
  sdp::MaxX   = MaxX;

  // Individual
  sdp::B0 = B0;
  sdp::w = w;
  sdp::xc = xc;
  sdp::expFor = expFor,
  sdp::expFly = expFly,
  sdp::b0 = b0;
  sdp::b1 = b1;
  sdp::b2 = b2;
  sdp::pred_a1 = pred_a1;
  sdp::pred_a2 = pred_a2;
  sdp::c = c;
  sdp::speed = speed;
  sdp::WindAssist = WindAssist;
  sdp::WindProb = WindProb;
  sdp::ZStdNorm = ZStdNorm;
  sdp::PStdNorm = PStdNorm;
  sdp::nTR_x = nTR_x;
  sdp::nTR_y = nTR_y;
  sdp::decError = decError;

  // Sites
  sdp::dist   = dist;
  sdp::bear   = bear;
  sdp::y_gain = y_gain;
  sdp::expend = y_expend;

}

// [[Rcpp::export]]
arma::vec simForaging(double f_intensity, int time, int site, int x, int inf)
{
  double gain, new_x = 0.00, mean, SD, Predat, Runif, pre_x;
  double LB, UB;
  int hit = 0, dead = false;
  int NStdNorm = sdp::ZStdNorm.size();
  arma::vec out = arma::zeros<arma::vec>(2);

  Runif = R::runif(0,1);
  LB = 0.0;
  UB = 0.0;
  for (int h = 0; h < NStdNorm; ++h)
  {
    UB += sdp::PStdNorm(h);
    if ((Runif > LB) && (Runif <= UB))
      hit = h;
    LB += sdp::PStdNorm(h);
  }
  mean = sdp::y_gain(site);
  SD =   1;

  double temp_gain = mean + sdp::ZStdNorm(hit) * SD;
  gain = temp_gain;
  if (gain < 0.0) gain = 0.0;

  pre_x = double(x);
  double expenditure = sdp::expend(site, time);

  new_x = pre_x + f_intensity * gain - expenditure;

  switch(Eval(new_x, 0.0, sdp::MaxX))
  {
  case 0: x = 0; break;
  case 1: x = Round(new_x); break;
  case 2: x = sdp::MaxX; break;
  };

  if ((f_intensity > 0.0) && (gain > 0.0) && (pre_x > 0.0))
  {
    double help1 = exp((sdp::pred_a1 + 1.0) * log(pre_x + f_intensity * gain));
    double help2 = exp((sdp::pred_a1 + 1.0) * log(pre_x));
    double help3 = (sdp::pred_a1 + 1.0) * f_intensity * gain;
    double help4 = exp(sdp::pred_a2 * log(f_intensity));
    Predat = sdp::b0(site) + sdp::b1(site) *
      ((help1 - help2)/help3) * sdp::b2(site) * help4;
  }
  else Predat = sdp::b0(site);

  if (Predat > 0.0)
  {
    if (R::runif(0, 1) < Predat) dead = true;
  }

  out(0) = Round(new_x);
  out(1) = dead;
  return out;
}

// [[Rcpp::export]]
arma::vec simFlying(int decision, int time, int site, int x, int inf)
{
  double nextx, total_D, distance, range;
  double Sqr_c, Sqr_ca, t, UB, LB, Runif;
  int hit = 0;
  int NWind = sdp::WindAssist.size();
  arma::vec out = arma::zeros<arma::vec>(2);

  int dest_site = decision;
  total_D = 0.0;
  total_D = sdp::dist(site, dest_site);

  Runif = R::runif(0,1);
  LB = 0.0;
  UB = 0.0;
  for (int h = 0; h < NWind; ++h)
  {
    UB += sdp::WindProb(h);
    if ((Runif > LB ) && (Runif <= UB))
      hit = h;
    LB += sdp::WindProb(h);
  }
  distance = total_D * (1.0 + sdp::WindAssist(hit));

  range =   (sdp::c) * (1.0 - (1.0/ (sqrt(1.0 + ( ((double)(x)) / (double)(sdp::MaxX))))));
  Sqr_c  = ((sdp::c) *(sdp::c));
  Sqr_ca = ((sdp::c) - (range - distance))*((sdp::c) - (range - distance));
  nextx  = (((Sqr_c/Sqr_ca) - 1.0) * (double)(sdp::MaxX));


  t =  time + (distance/ sdp::speed);

  out(0) = Round(t);
  out(1) = Round(nextx);
  return out;
}


