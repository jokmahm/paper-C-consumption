#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()(){
  // data----------------------------------------------------------
  DATA_INTEGER(ny);
  DATA_SCALAR(T_0);
  DATA_SCALAR(T_l);
  
  DATA_VECTOR(P);
  DATA_VECTOR(B);
  DATA_VECTOR(Temp);
  DATA_VECTOR(f);
  
  DATA_VECTOR(predP);
  DATA_VECTOR(predB);
  DATA_VECTOR(predTemp);
  DATA_INTEGER(logbhTRUE);
  
  // parameters----------------------------------------------------
  PARAMETER(x); 
  PARAMETER(logb_alpha); 
  PARAMETER(logb_h);
  PARAMETER(BH);
  PARAMETER(logSigma);
  
  
  //Type x = exp(logx);
  Type b_alpha = exp(logb_alpha);
  Type b_h =0;
  if(logbhTRUE == 1){
    b_h = exp(logb_h);
  }else{
    b_h = BH;
  }
  Type sigma = exp(logSigma); 
  // Type xtrans = Type(2.0) * exp(x)/(Type(1.0) + exp(x))-Type(1.0);
  
  // Define needed variables---------------------------------------
  Type nll=0;
  
  //Define some needed vectors h and alpha ------------------------
  vector<Type> h(ny);
  vector<Type> alpha(ny);
  
  //----------------------------------------------------------------
  vector<Type> C(ny);
  // Compute Regression---------------------------------------------
  for (int i=0;i<ny;i++) {
    alpha(i) = b_alpha * (Temp(i) - T_0) * pow((T_l - Temp(i)), Type(0.0));
    h(i)     = b_h; //* (Temp(i) - T_0) * pow((T_l-Temp(i)), Type(0.0));
    C(i) = (alpha(i)*pow(P(i),-x) * B(i)/(Type(1)+alpha(i)*h(i)*pow(P(i),-x)*B(i)));
    nll -= -log(f(i)) - log(sigma) - Type(0.5) * log(2.0 * M_PI) - pow((log(f(i)) - log(C(i))), 2) / (2.0 * sigma * sigma);
  }
  
  // Predict the remaining years
  vector<Type> predC(predB.size());
  for(int i =0; i < predB.size(); i++){
    Type altemp = b_alpha * (predTemp(i) - T_0) * pow((T_l - predTemp(i)), Type(0.0));
    Type htemp  = b_h; //* (predTemp(i) - T_0) * pow((T_l-predTemp(i)), Type(0.0));
    predC(i) = (altemp * pow(predP(i),-x) * predB(i) / (Type(1)+altemp *htemp * pow(predP(i),-x)*predB(i)));
  }
  
  
  
  
  //nll = sqrt(nll);
  ADREPORT(x);
  ADREPORT(b_alpha);
  ADREPORT(b_h);
  ADREPORT(sigma);
  ADREPORT(alpha);
  ADREPORT(h);
  ADREPORT(C);
  ADREPORT(predC);
  
  //ADREPORT(alpha);
  //ADREPORT(x);
  //ADREPORT(theta_0);
  //ADREPORT(theta_1);
  //ADREPORT(theta_2);
  return(nll);
}
