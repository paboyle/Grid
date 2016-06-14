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

#include <Hadrons/Global.hpp>
#include <Hadrons/Graph.hpp>

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
        for (auto &p: s.population_)
        {
            out << p.second << ": " << p.first << std::endl;
        }
        
        return out;
    }
private:
    // randomly initialize population
    void initPopulation(void);
    // genetic operators
    const std::vector<T> & selection(void);
    void                   crossover(const std::vector<T> &c1,
                                     const std::vector<T> &c2);
    void                   mutation(std::vector<T> &c);
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
    std::uniform_real_distribution<double> dis(0., 1.);
    
    // random initialization of the population if necessary
    if (population_.size() != par_.popSize)
    {
        initPopulation();
    }

    // mating
    for (unsigned int i = 0; i < par_.popSize/2; ++i)
    {
        auto &p1 = selection(), &p2 = selection();
        crossover(p1, p2);
    }
    
    // random mutations
    auto buf = population_;
    population_.clear();
    for (auto &c: buf)
    {
        if (dis(gen_) < par_.mutationRate)
        {
            mutation(c.second);
        }
        population_.emplace(func_(c.second), c.second);
    }
    
    // grim reaper
    auto it = population_.begin();
    
    std::advance(it, par_.popSize);
    population_.erase(it, population_.end());
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
const std::vector<T> & GeneticScheduler<T>::selection(void)
{
    std::vector<double> prob;
    
    for (auto &c: population_)
    {
        prob.push_back(1./c.first);
    }
    std::discrete_distribution<unsigned int> dis(prob.begin(), prob.end());
    auto rIt = population_.begin();
    std::advance(rIt, dis(gen_));
    
    return rIt->second;
}

template <typename T>
void GeneticScheduler<T>::crossover(const std::vector<T> &p1,
                                    const std::vector<T> &p2)
{
    std::uniform_int_distribution<unsigned int> dis(0, p1.size() - 1);
    unsigned int                                cut = dis(gen_);
    std::vector<T>                              c1, c2, buf;
    
    auto cross = [&buf, cut](std::vector<T> &c, const std::vector<T> &p1,
                             const std::vector<T> &p2)
    {
        buf = p2;
        for (unsigned int i = 0; i < cut; ++i)
        {
            c.push_back(p1[i]);
            buf.erase(std::find(buf.begin(), buf.end(), p1[i]));
        }
        for (unsigned int i = 0; i < buf.size(); ++i)
        {
            c.push_back(buf[i]);
        }
    };
    
    cross(c1, p1, p2);
    cross(c2, p2, p1);
    population_.emplace(func_(c1), c1);
    population_.emplace(func_(c2), c2);
}

template <typename T>
void GeneticScheduler<T>::mutation(std::vector<T> &c)
{
    std::uniform_int_distribution<unsigned int> dis(0, c.size() - 1);
    unsigned int                                cut = dis(gen_);
    Graph<T>                                    g = graph_;
    std::vector<T>                              buf;
    
    for (unsigned int i = cut; i < c.size(); ++i)
    {
        g.removeVertex(c[i]);
    }
    if (g.size() > 0)
    {
        buf = g.topoSort(gen_);
    }
    for (unsigned int i = cut; i < c.size(); ++i)
    {
        buf.push_back(c[i]);
    }
    c = buf;
}

END_HADRONS_NAMESPACE

#endif // Hadrons_GeneticScheduler_hpp_
