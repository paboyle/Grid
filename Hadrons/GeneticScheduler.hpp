/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/GeneticScheduler.hpp

Copyright (C) 2015-2019

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_GeneticScheduler_hpp_
#define Hadrons_GeneticScheduler_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Graph.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Scheduler based on a genetic algorithm                   *
 ******************************************************************************/
template <typename V, typename T>
class GeneticScheduler
{
public:
    typedef std::vector<T>                 Gene;
    typedef std::pair<Gene *, Gene *>      GenePair;
    typedef std::function<V(const Gene &)> ObjFunc;
    struct Parameters
    {
        double       mutationRate;
        unsigned int popSize, seed;
    };
public:
    // constructor
    GeneticScheduler(Graph<T> &graph, const ObjFunc &func,
                     const Parameters &par);
    // destructor
    virtual ~GeneticScheduler(void) = default;
    // access
    const Gene & getMinSchedule(void);
    V            getMinValue(void);
    // reset population
    void initPopulation(void);
    // breed a new generation
    void nextGeneration(void);
    // heuristic benchmarks
    void benchmarkCrossover(const unsigned int nIt);
    // print population
    friend std::ostream & operator<<(std::ostream &out,
                                     const GeneticScheduler<V, T> &s)
    {
        out << "[";
        for (auto &p: s.population_)
        {
            out << p.first << ", ";
        }
        out << "\b\b]";
        
        return out;
    }
private:
    void doCrossover(void);
    void doMutation(void);
    // genetic operators
    GenePair selectPair(void);
    void     crossover(Gene &c1, Gene &c2, const Gene &p1, const Gene &p2);
    void     mutation(Gene &m, const Gene &c);
    
private:
    Graph<T>               &graph_;
    const ObjFunc          &func_;
    const Parameters       par_;
    std::multimap<V, Gene> population_;
    std::mt19937           gen_;
};

/******************************************************************************
 *                       template implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename V, typename T>
GeneticScheduler<V, T>::GeneticScheduler(Graph<T> &graph, const ObjFunc &func,
                                      const Parameters &par)
: graph_(graph)
, func_(func)
, par_(par)
{
    gen_.seed(par_.seed);
}

// access //////////////////////////////////////////////////////////////////////
template <typename V, typename T>
const typename GeneticScheduler<V, T>::Gene &
GeneticScheduler<V, T>::getMinSchedule(void)
{
    return population_.begin()->second;
}

template <typename V, typename T>
V GeneticScheduler<V, T>::getMinValue(void)
{
    return population_.begin()->first;
}

// breed a new generation //////////////////////////////////////////////////////
template <typename V, typename T>
void GeneticScheduler<V, T>::nextGeneration(void)
{
    // random initialization of the population if necessary
    if (population_.size() != par_.popSize)
    {
        initPopulation();
    }
    //LOG(Debug) << "Starting population:\n" << *this << std::endl;
    
    // random mutations
    //PARALLEL_FOR_LOOP
    for (unsigned int i = 0; i < par_.popSize; ++i)
    {
        doMutation();
    }
    //LOG(Debug) << "After mutations:\n" << *this << std::endl;
    
    // mating
    //PARALLEL_FOR_LOOP
    for (unsigned int i = 0; i < par_.popSize/2; ++i)
    {
        doCrossover();
    }
    //LOG(Debug) << "After mating:\n" << *this << std::endl;
    
    // grim reaper
    auto it = population_.begin();
    
    std::advance(it, par_.popSize);
    population_.erase(it, population_.end());
    //LOG(Debug) << "After grim reaper:\n" << *this << std::endl;
}

// evolution steps /////////////////////////////////////////////////////////////
template <typename V, typename T>
void GeneticScheduler<V, T>::initPopulation(void)
{
    population_.clear();
    for (unsigned int i = 0; i < par_.popSize; ++i)
    {
        auto p = graph_.topoSort(gen_);
        
        population_.insert(std::make_pair(func_(p), p));
    }
}

template <typename V, typename T>
void GeneticScheduler<V, T>::doCrossover(void)
{
    auto p = selectPair();
    Gene &p1 = *(p.first), &p2 = *(p.second);
    Gene c1, c2;
    
    crossover(c1, c2, p1, p2);
    PARALLEL_CRITICAL
    {
        population_.insert(std::make_pair(func_(c1), c1));
        population_.insert(std::make_pair(func_(c2), c2));
    }
}

template <typename V, typename T>
void GeneticScheduler<V, T>::doMutation(void)
{
    std::uniform_real_distribution<double>      mdis(0., 1.);
    std::uniform_int_distribution<unsigned int> pdis(0, population_.size() - 1);
    
    if (mdis(gen_) < par_.mutationRate)
    {
        Gene m;
        auto it = population_.begin();
        
        std::advance(it, pdis(gen_));
        mutation(m, it->second);
        PARALLEL_CRITICAL
        {
            population_.insert(std::make_pair(func_(m), m));
        }
    }
}

// genetic operators ///////////////////////////////////////////////////////////
template <typename V, typename T>
typename GeneticScheduler<V, T>::GenePair GeneticScheduler<V, T>::selectPair(void)
{
    std::vector<double> prob;
    unsigned int        ind;
    Gene                *p1, *p2;
    const double        max = population_.rbegin()->first;
    

    for (auto &c: population_)
    {
        prob.push_back(std::exp((c.first-1.)/max));
    }        
    std::discrete_distribution<unsigned int> dis1(prob.begin(), prob.end());
    auto rIt = population_.begin();
    ind = dis1(gen_);
    std::advance(rIt, ind);
    p1 = &(rIt->second);
    prob[ind] = 0.;
    std::discrete_distribution<unsigned int> dis2(prob.begin(), prob.end());
    rIt = population_.begin();
    std::advance(rIt, dis2(gen_));
    p2 = &(rIt->second);
    
    return std::make_pair(p1, p2);
}

template <typename V, typename T>
void GeneticScheduler<V, T>::crossover(Gene &c1, Gene &c2, const Gene &p1,
                                    const Gene &p2)
{
    Gene                                        buf;
    std::uniform_int_distribution<unsigned int> dis(0, p1.size() - 1);
    unsigned int                                cut = dis(gen_);
    
    c1.clear();
    buf = p2;
    for (unsigned int i = 0; i < cut; ++i)
    {
        c1.push_back(p1[i]);
        buf.erase(std::find(buf.begin(), buf.end(), p1[i]));
    }
    for (unsigned int i = 0; i < buf.size(); ++i)
    {
        c1.push_back(buf[i]);
    }
    c2.clear();
    buf = p2;
    for (unsigned int i = cut; i < p1.size(); ++i)
    {
        buf.erase(std::find(buf.begin(), buf.end(), p1[i]));
    }
    for (unsigned int i = 0; i < buf.size(); ++i)
    {
        c2.push_back(buf[i]);
    }
    for (unsigned int i = cut; i < p1.size(); ++i)
    {
        c2.push_back(p1[i]);
    }
}

template <typename V, typename T>
void GeneticScheduler<V, T>::mutation(Gene &m, const Gene &c)
{
    Gene                                        buf;
    std::uniform_int_distribution<unsigned int> dis(0, c.size() - 1);
    unsigned int                                cut = dis(gen_);
    Graph<T>                                    g1 = graph_, g2 = graph_;
    
    for (unsigned int i = 0; i < cut; ++i)
    {
        g1.removeVertex(c[i]);
    }
    for (unsigned int i = cut; i < c.size(); ++i)
    {
        g2.removeVertex(c[i]);
    }
    if (g1.size() > 0)
    {
        buf = g1.topoSort(gen_);
    }
    if (g2.size() > 0)
    {
        m = g2.topoSort(gen_);
    }
    for (unsigned int i = cut; i < c.size(); ++i)
    {
        m.push_back(buf[i - cut]);
    }
}

template <typename V, typename T>
void GeneticScheduler<V, T>::benchmarkCrossover(const unsigned int nIt)
{
    Gene   p1, p2, c1, c2;
    double neg = 0., eq = 0., pos = 0., total;
    int    improvement;
    
    LOG(Message) << "Benchmarking crossover..." << std::endl;
    for (unsigned int i = 0; i < nIt; ++i)
    {
        p1 = graph_.topoSort(gen_);
        p2 = graph_.topoSort(gen_);
        crossover(c1, c2, p1, p2);
        improvement = (func_(c1) + func_(c2) - func_(p1) - func_(p2))/2;
        if (improvement < 0) neg++; else if (improvement == 0) eq++; else pos++;
    }
    total = neg + eq + pos;
    LOG(Message) << "  -: " << neg/total << "  =: " << eq/total
                 << "  +: " << pos/total << std::endl;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_GeneticScheduler_hpp_
