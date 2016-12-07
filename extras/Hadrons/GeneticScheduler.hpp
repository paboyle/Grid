/*******************************************************************************
Grid physics library, www.github.com/paboyle/Grid 

Source file: programs/Hadrons/GeneticScheduler.hpp

Copyright (C) 2016

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

See the full license in the file "LICENSE" in the top level distribution 
directory.
*******************************************************************************/

#ifndef Hadrons_GeneticScheduler_hpp_
#define Hadrons_GeneticScheduler_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Graph.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                   Scheduler based on a genetic algorithm                   *
 ******************************************************************************/
template <typename T>
class GeneticScheduler
{
public:
    typedef std::function<int(const std::vector<T> &)> ObjFunc;
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
    const std::vector<T> & getMinSchedule(void);
    int                    getMinValue(void);
    // breed a new generation
    void nextGeneration(void);
    // print population
    friend std::ostream & operator<<(std::ostream &out,
                                     const GeneticScheduler<T> &s)
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
    // randomly initialize population
    void initPopulation(void);
    // genetic operators
    std::vector<T> *                              select1(void);
    std::pair<std::vector<T> *, std::vector<T> *> select2(void);
    void                                          crossover(void);
    void                                          mutation(void);
private:
    Graph<T>                           &graph_;
    const ObjFunc                      &func_;
    const Parameters                   par_;
    std::multimap<int, std::vector<T>> population_;
    std::mt19937                       gen_;
};

/******************************************************************************
 *                       template implementation                              *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
GeneticScheduler<T>::GeneticScheduler(Graph<T> &graph, const ObjFunc &func,
                                      const Parameters &par)
: graph_(graph)
, func_(func)
, par_(par)
{
    gen_.seed(par_.seed);
}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
const std::vector<T> & GeneticScheduler<T>::getMinSchedule(void)
{
    return population_.begin()->second;
}

template <typename T>
int GeneticScheduler<T>::getMinValue(void)
{
    return population_.begin()->first;
}

// breed a new generation //////////////////////////////////////////////////////
template <typename T>
void GeneticScheduler<T>::nextGeneration(void)
{
    // random initialization of the population if necessary
    if (population_.size() != par_.popSize)
    {
        initPopulation();
    }
    LOG(Debug) << "Starting population:\n" << *this << std::endl;
    
    // random mutations
    PARALLEL_FOR_LOOP
    for (unsigned int i = 0; i < par_.popSize; ++i)
    {
        mutation();
    }
    LOG(Debug) << "After mutations:\n" << *this << std::endl;
    
    // mating
    PARALLEL_FOR_LOOP
    for (unsigned int i = 0; i < par_.popSize/2; ++i)
    {
        crossover();
    }
    LOG(Debug) << "After mating:\n" << *this << std::endl;
    
    // grim reaper
    auto it = population_.begin();
    
    std::advance(it, par_.popSize);
    population_.erase(it, population_.end());
    LOG(Debug) << "After grim reaper:\n" << *this << std::endl;
}

// randomly initialize population //////////////////////////////////////////////
template <typename T>
void GeneticScheduler<T>::initPopulation(void)
{
    population_.clear();
    for (unsigned int i = 0; i < par_.popSize; ++i)
    {
        auto p = graph_.topoSort(gen_);
        
        population_.emplace(func_(p), p);
    }
}

// genetic operators ///////////////////////////////////////////////////////////
template <typename T>
std::vector<T> * GeneticScheduler<T>::select1(void)
{
    std::uniform_int_distribution<unsigned int> pdis(0, population_.size() - 1);
    
    auto it = population_.begin();
    std::advance(it, pdis(gen_));
    
    return &(it->second);
}

template <typename T>
std::pair<std::vector<T> *, std::vector<T> *> GeneticScheduler<T>::select2(void)
{
    std::vector<double> prob;
    unsigned int        ind;
    std::vector<T>      *p1, *p2;
    
    for (auto &c: population_)
    {
        prob.push_back(1./c.first);
    }
    do
    {
        double probCpy;
        
        std::discrete_distribution<unsigned int> dis1(prob.begin(), prob.end());
        auto rIt = population_.begin();
        ind = dis1(gen_);
        std::advance(rIt, ind);
        p1 = &(rIt->second);
        probCpy   = prob[ind];
        prob[ind] = 0.;
        std::discrete_distribution<unsigned int> dis2(prob.begin(), prob.end());
        rIt = population_.begin();
        std::advance(rIt, dis2(gen_));
        p2 = &(rIt->second);
        prob[ind] = probCpy;
    } while (p1 == p2);
    
    return std::make_pair(p1, p2);
}

template <typename T>
void GeneticScheduler<T>::crossover(void)
{
    auto                                        p = select2();
    auto                                        &p1 = *(p.first),
                                                &p2 = *(p.second);
    std::uniform_int_distribution<unsigned int> dis2(0, p1.size() - 1);
    unsigned int                                cut = dis2(gen_);
    std::vector<T>                              c1, c2, buf;
    
    auto cross = [&buf, cut](std::vector<T> &c, const std::vector<T> &p1,
                             const std::vector<T> &p2)
    {
        buf = p1;
        for (unsigned int i = 0; i < cut; ++i)
        {
            c.push_back(p2[i]);
            buf.erase(std::find(buf.begin(), buf.end(), p2[i]));
        }
        for (unsigned int i = 0; i < buf.size(); ++i)
        {
            c.push_back(buf[i]);
        }
    };
    
    cross(c1, p1, p2);
    cross(c2, p2, p1);
    PARALLEL_CRITICAL
    {
        population_.emplace(func_(c1), c1);
        population_.emplace(func_(c2), c2);
    }
}

template <typename T>
void GeneticScheduler<T>::mutation(void)
{
    std::uniform_real_distribution<double> mdis(0., 1.);
    
    if (mdis(gen_) < par_.mutationRate)
    {
        auto                                        &c = *select1();
        std::uniform_int_distribution<unsigned int> cdis(0, c.size() - 1);
        unsigned int                                cut = cdis(gen_);
        std::vector<T>                              buf1, buf2;
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
            buf1 = g1.topoSort(gen_);
        }
        if (g2.size() > 0)
        {
            buf2 = g2.topoSort(gen_);
        }
        for (unsigned int i = cut; i < c.size(); ++i)
        {
            buf2.push_back(buf1[i - cut]);
        }
        PARALLEL_CRITICAL
        {
            population_.emplace(func_(buf2), buf2);
        }
    }
}

END_HADRONS_NAMESPACE

#endif // Hadrons_GeneticScheduler_hpp_
