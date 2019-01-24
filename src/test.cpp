//============================================================================
// Name        : test.cpp
// Author      : f
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fst/fst-decl.h>
#include <fst/dfs-visit.h>
#include <fst/visit.h>
#include <fstext/fstext-utils.h>
#include <fstext/table-matcher.h>
#include <base/timer.h>

using namespace std;
using namespace fst;
using namespace kaldi;

ComposeFst<StdArc>* OTFComposeFst(
  const StdFst &ifst1, const StdFst &ifst2,
  const CacheOptions& cache_opts = CacheOptions()) {

  typedef LookAheadMatcher< StdFst > M;
  typedef AltSequenceComposeFilter<M> SF;
  typedef LookAheadComposeFilter<SF, M>  LF;
  typedef PushWeightsComposeFilter<LF, M> WF;
  typedef PushLabelsComposeFilter<WF, M> ComposeFilter;
  typedef M FstMatcher;

  ComposeFstOptions<StdArc, FstMatcher, ComposeFilter> opts(cache_opts);

  return new ComposeFst<StdArc>(ifst1, ifst2, opts);
//    return new ComposeFst<StdArc>(ifst1, ifst2);
}

ComposeFst<StdArc> TableComposeFst(const StdFst &ifst1, const StdFst &ifst2, CacheOptions &cache_opts) {
  typedef SortedMatcher<StdFst> SM;
  typedef TableMatcher<StdFst> TM;
  typedef ArcLookAheadMatcher<SM> LA_SM;
  typedef SequenceComposeFilter<TM, LA_SM> SCF;
  typedef LookAheadComposeFilter<SCF, TM, LA_SM, MATCH_INPUT> LCF;
  typedef PushWeightsComposeFilter<LCF, TM, LA_SM, MATCH_INPUT> PWCF;
  typedef PushLabelsComposeFilter<PWCF, TM, LA_SM, MATCH_INPUT> PWLCF;
  typedef GenericComposeStateTable<StdArc,PWLCF::FilterState> ST;
  TM* lam1 = new TM(ifst1, MATCH_OUTPUT);
  LA_SM* lam2 = new LA_SM(ifst2, MATCH_INPUT);
  PWLCF* laf = new PWLCF(ifst1, ifst2, lam1, lam2);
  ST* st = new ST(ifst1, ifst2);
  ComposeFstImplOptions<TM, LA_SM, PWLCF> opts(cache_opts, lam1, lam2, laf, st);
  return ComposeFst<StdArc>(ifst1, ifst2, opts);
}

typedef VectorFst<LogArc> LogVectorFst;

int pathCount(const StdFst &sfst) {
	// if the fst has loops, then the language or relation
	// encoded by the FST is infinite, just return -1
	if (sfst.Properties(kCyclic, true)) {
		return -1;
	}
	// copy the FST input to a new unweighted network,
	// over the Log Semiring
	LogVectorFst lfst;

	// for each original StdArc, if the weight is _not_
	// TropicalWeight::Zero() then change it to LogWeight::One(),
	// else (when the weight is TropicalWeight::Zero()) change it to
	// LogWeight::Zero().  This is handily done with
	// RmWeightMapper(), which is in the library.

	Map(sfst, &lfst, RmWeightMapper<StdArc, LogArc>());

	vector<LogArc::Weight> distance;
	ShortestDistance(lfst, &distance, true);

	// from the distance vector, get the weight for the start state
	LogArc::Weight w = distance[lfst.Start()] ;

	// w.Value() is the -log of the number of paths
	// cout << "Weight from pathCount(): " << w << endl ;
	// so make w positive and get the exp()
	int paths = exp((double)(-1.0 * w.Value()));
	// cout << "Paths from pathCount(): " << paths << endl;
	return paths ;
}

typedef StdArc::StateId StateId;

void dfs(StdFst &fst, StateId s, vector<bool> &visited) {
	visited[s] = true;
	for (ArcIterator<StdFst> aiter(fst, s); !aiter.Done(); aiter.Next()) {
		const StdArc &arc = aiter.Value();
		if (!visited[arc.nextstate]) dfs(fst, arc.nextstate, visited);
	}
}

typedef SortedMatcher<StdFst> SM;
typedef TableMatcher<StdFst> TM;
typedef ArcLookAheadMatcher<SM> LA_SM;
typedef SequenceComposeFilter<TM, LA_SM> SCF;
typedef LookAheadComposeFilter<SCF, TM, LA_SM, MATCH_INPUT> LCF;
typedef PushWeightsComposeFilter<LCF, TM, LA_SM, MATCH_INPUT> PWCF;
typedef PushLabelsComposeFilter<PWCF, TM, LA_SM, MATCH_INPUT> PWLCF;
typedef GenericComposeStateTable<StdArc,PWLCF::FilterState> ST;
typedef CacheImplOptions<DefaultCacheStore<StdArc>> CACHE_OPTIONS;

ComposeFst<StdArc> OTFCompose(StdFst &ifst1, StdFst &ifst2, CACHE_OPTIONS &cache_opts)
{
	TM* lam1 = new TM(ifst1, MATCH_OUTPUT);
	LA_SM* lam2 = new LA_SM(ifst2, MATCH_INPUT);
	PWLCF* laf = new PWLCF(ifst1, ifst2, lam1, lam2);
	ST* st = new ST(ifst1, ifst2);
	ComposeFstImplOptions<TM, LA_SM, PWLCF> opts(cache_opts, lam1, lam2, laf, st);
	return ComposeFst<StdArc>(ifst1, ifst2, opts);
}

int main() {

	StdFst *left_fst = Fst<StdArc>::Read("otf_gru/graph_V1_lm1_nnet3/HCLG_mod.fst");
	StdFst *right_fst = Fst<StdArc>::Read("otf_gru/lang_V+Add_lm3/replace.fst");

//	ComposeFst<StdArc> compose_fst = TableComposeFst(*left_fst, *right_fst);
	//ComposeFst<StdArc> compose_fst = TableComposeFst(*left_fst, *right_fst);

	DefaultCacheStore<StdArc> store =
			DefaultCacheStore<StdArc>(CacheOptions(false,2000000));
	CACHE_OPTIONS cache_opts = CACHE_OPTIONS(false, 2000000, &store);

	cache_opts.own_store = false;

	//ComposeFst<StdArc> compose_fst =  TableComposeFst(*left_fst, *right_fst, cache_opts);

	for (int i=0; i<2; i++)
	{
		ComposeFst<StdArc> compose_fst = OTFCompose(*left_fst, *right_fst, cache_opts);
		Timer timer;
		int state_count = 0;
		for (StateIterator<StdFst> siter(compose_fst); !siter.Done(); siter.Next())
			cout << state_count++ << endl;
		vector<bool> visited(state_count);
		for (auto e : visited) e = false;
		dfs(compose_fst, compose_fst.Start(), visited);
		cout << "Time elapsed " << timer.Elapsed() << endl;
	}

//	if (compose_fst.Write("out.fst")) cout << "good" << endl;
//	else cout << "bad" << endl;

	//int count = pathCount(compose_fst);
	//cout << count << endl;

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
