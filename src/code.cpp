#include <Rcpp.h>
//# © 2023 Rutgers, The State University of New Jersey. All rights reserved.


using namespace Rcpp;


inline int randWrapper(const int n) { return floor(unif_rand()*n); }

NumericVector randomShuffle2(NumericVector a) {
  // clone a into b to leave a alone
  NumericVector b =clone(a);
  int n = b.size();
  int j;
  
  // Fisher-Yates Shuffle Algorithm
  for (int i = 0; i < n - 1; i++) {
    j = i + randWrapper(n - i);
    std::swap(b[i], b[j]);
  }
  return b;
}

NumericVector sort_custom(NumericVector x) {
  IntegerVector idx = seq_along(x) - 1;
  std::sort(idx.begin(), idx.end(), [&](int i, int j){return x[i] > x[j];});
  return x[idx];
}


//' Single sample GSEA
//' 
//' This function returns a List basically a list with the
//' following elements:
//' @return \item{ES}{Enrichiment Score of the pathway w.r.t to gene}
//' @return \item{NES}{Normailized Enrinchment Score}
//' @return \item{p.value}{P value of the estimate}
//' @references 
// [[Rcpp::export]]
List ssgsea(NumericVector ref,CharacterVector gs, int np=100){
  NumericVector reflist  = sort_custom(ref);
  
  CharacterVector reflist_names = reflist.names();
  
  gs= intersect(reflist_names,gs);
  double es= 0.0;
  double nes= 0.0;
  double pval= 0.0;
  double p_val_sum =0.0;
  double p_val_denm =0;
  
  
  // null value checks 
  if(!is_true(any(is_na(gs))) and !is_true(any(is_na(reflist))) and gs.size()>=1){
    NumericVector isgs(reflist.size());
    
    LogicalVector cond = Rcpp::in(reflist_names, gs);
    for(unsigned int i = 0; i < cond.size(); ++i){
      if(cond[i]) {
        isgs[i] = 1;
      } 
    }
    
    
    int n = reflist.size();
    NumericVector score_hit(n);
    double sum = 0;
    
    for (int i = 0; i < n; i++) {
      double product = reflist[i] * isgs[i];
      sum += std::abs(product);
      score_hit[i] = sum;
    }
    
    
    score_hit = score_hit/tail(score_hit, 1)[0];
    NumericVector score_miss = cumsum(1-isgs);
    score_miss = score_miss/tail(score_miss, 1)[0];
    NumericVector score_all = score_hit - score_miss;
    es = max(score_all) + min(score_all);
    
    
    double total_bg_es= 0.0;
    int count_bg_es =0;

    for(int i = 0; i < np; ++i) {
      NumericVector bg_isgs = randomShuffle2(isgs);
      NumericVector bg_hit = cumsum(abs(reflist*bg_isgs));
      bg_hit = bg_hit/tail(bg_hit, 1)[0];
      NumericVector bg_miss = cumsum(1-bg_isgs);
      bg_miss = bg_miss/tail(bg_miss, 1)[0];
      NumericVector bg_all = bg_hit - bg_miss;
      double bges = max(bg_all) + min(bg_all);
      
      if(es<0){
        if(bges<0){
          p_val_denm+=1;
        }
        if(bges<=es){
          p_val_sum+=1;
          
        }
      }else{
        if(bges>0){
          p_val_denm+=1;
        }
        if(bges>=es){
          p_val_sum+=1;
        }
      }
      
      
      if(es<0 && bges <0){
        total_bg_es+=bges;
        count_bg_es+=1;
      }
      if(es>0 && bges >0){
        total_bg_es+=bges;
        count_bg_es+=1;
      }
      
    }
    
    nes = 0;
    pval = 1;
    if(abs(total_bg_es/count_bg_es)>0)
    {
      nes = es/abs(total_bg_es/count_bg_es);
    }
    
  if (p_val_denm > 0)
  {
    pval = p_val_sum/p_val_denm;
    }
  
    
  }

  if (pval ==1.0){
    pval = 1/np;
  }
  
  List L = List::create(Named("ES") = es , _["NES"] = nes,_["p.value"] = pval);
  return L;
}
