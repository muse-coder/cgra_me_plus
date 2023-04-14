/*******************************************************************************
 * CGRA-ME Software End-User License Agreement
 *
 * The software programs comprising "CGRA-ME" and the documentation provided
 * with them are copyright by its authors S. Chin, K. Niu, N. Sakamoto, J. Zhao,
 * A. Rui, S. Yin, A. Mertens, J. Anderson, and the University of Toronto. Users
 * agree to not redistribute the software, in source or binary form, to other
 * persons or other institutions. Users may modify or use the source code for
 * other non-commercial, not-for-profit research endeavours, provided that all
 * copyright attribution on the source code is retained, and the original or
 * modified source code is not redistributed, in whole or in part, or included
 * in or with any commercial product, except by written agreement with the
 * authors, and full and complete attribution for use of the code is given in
 * any resulting publications.
 *
 * Only non-commercial, not-for-profit use of this software is permitted. No
 * part of this software may be incorporated into a commercial product without
 * the written consent of the authors. The software may not be used for the
 * design of a commercial electronic product without the written consent of the
 * authors. The use of this software to assist in the development of new
 * commercial CGRA architectures or commercial soft processor architectures is
 * also prohibited without the written consent of the authors.
 *
 * This software is provided "as is" with no warranties or guarantees of
 * support.
 *
 * This Agreement shall be governed by the laws of Province of Ontario, Canada.
 *
 * Please contact Prof. Anderson if you are interested in commercial use of the
 * CGRA-ME framework.
 ******************************************************************************/

#ifndef ___ANNEALMAPPER_H__
#define ___ANNEALMAPPER_H__

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <algorithm>
#include <unordered_set>
#include <CGRA/CGRA.h>
#include <CGRA/OpGraph.h>
#include <CGRA/Mapper.h>
#include <CGRA/Mapping.h>

using namespace std;





class AnnealMapper : public Mapper
{
public:
    // Custom setup with variable args
    AnnealMapper(std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> &args);
    // Default setup
    AnnealMapper(std::shared_ptr<CGRA> cgra);
    // Custom setup
    AnnealMapper(std::shared_ptr<CGRA> cgra, int timelimit, int rand_seed, float initial_penalty, float penalty_factor, float const_temp_factor, int swap_factor, float cold_accept_rate);
    ~AnnealMapper() = default;

    Mapping mapOpGraph(std::shared_ptr<OpGraph> opgraph, int II, std::string arch_model_name) override;

    static void topological_sort_visit(std::vector<OpGraphOp *> *L, std::map<OpGraphOp *, bool> *mark, OpGraphOp *n)
    {
        if (!(*mark)[n])
        {
            if (n->opcode != OPGRAPH_OP_OUTPUT)
            {
                for (auto &m : n->output->output)
                {
                    topological_sort_visit(L, mark, m);
                }
            }

            (*mark)[n] = true;
            L->insert(L->begin(), n); // push n to head
        }
    }

    static OpGraphOp *topological_sort_check_marked(std::vector<OpGraphOp *> nodes, std::map<OpGraphOp *, bool> *mark)
    {
        for (auto &n : nodes)
        {
            if (!(*mark)[n])
                return n;
        }
        return NULL;
    }

    // This topological sort algorithm uses a depth first search
    // NB: The graph MUST be acyclic, else explosions may happen...
    void topological_sort(std::shared_ptr<OpGraph> g)
    {
        std::vector<OpGraphOp *> L;

        // // Unmark all nodes
        // std::map<OpGraphOp *, bool> mark;
        // for (auto &n : g->op_nodes)
        // {
        //     mark[n] = false;
        // }

        // OpGraphOp *n;
        // while ((n = topological_sort_check_marked(g->op_nodes, &mark)))
        // {
        //     topological_sort_visit(&L, &mark, n);
        // }
        std::set<OpGraphOp *> waiting_set;
         for (auto n : g->op_nodes)
        {
            waiting_set.insert(n);
        }
        int iteration = 0;
        while(iteration<1000){
            for(auto n: waiting_set){
                bool mapped  = true;
                for(auto input_val: n->input){
                    auto input_node = input_val->input;
                    if(mapping[input_node].size() == 0){
                        mapped = false;
                    }
                }
                if(mapped){
                    L.push_back(n);
                    waiting_set.erase(n);
                    break;
                }
            }

            iteration ++;

        }

        for(auto node: waiting_set){
            L.push_back(node);
        }

        assert(L.size() ==  g->op_nodes.size());
        g->op_nodes = L;
    }
    /*
bool test_compare(float penalty_factor, MRRGNode* a, MRRGNode* b)
{
    return getCost(a) < getCost(b);
}

bool compare_op_cost(float penalty_factor, OpGraphOp* a, OpGraphOp* b)
{
    return a->getCost(penalty_factor) > b->getCost(penalty_factor);
}
*/
    float updateTemperature(float t, float acceptance_rate)
    {
        if (acceptance_rate > 0.96)
        {
            return t * 0.5;
        }
        else if (acceptance_rate > 0.8)
        {
            return t * 0.9;
        }
        else if (acceptance_rate > 0.15)
        {
            return t * 0.98;
        }
        else
        {
            return t * 0.95;
        }
    }

    std::vector<std::pair<MRRGNode *, int>> candidate_fu_intersect(std::vector<std::pair<MRRGNode *, int>> v1, std::vector<std::pair<MRRGNode *, int>> v2)
    {
        std::vector<std::pair<MRRGNode *, int>> result;

        for (auto &p : v1)
        {
            for (auto &n : v2)
            {
                if (p.first == n.first)
                    result.push_back({p.first, p.second + n.second});
            }
        }

        return result;
    }

    // generate a random FU, zhaoying changes it to non-occupied
    MRRGNode *getRandomFU(MRRG *mrrg, OpGraphOp *op)
    {
        int min_start = 0;
        std::vector<MRRGNode *> candidates;
        for (auto &fu : mrrg->function_nodes)
        {   
            mpos mmmm{fu->x_, fu->y_, fu->cycle};
            bool exist = false;
            for(auto m :  mapped_pos){
                if(m.x == mmmm.x && mmmm.y == m.y  && mmmm.t == m.t){
                    exist = true;
                    break;
                }
            }
            if (fu->canMapOp(op)&& occupancy[fu] == 0)
            {
                candidates.push_back(fu);
            }
        }

        if (candidates.size() < 1)
        {
            cout << "Could not place: " << *op << endl;
            assert(candidates.size() > 0);
        }

        return candidates[rand() % candidates.size()];
    }

    template <typename type>
    std::vector<type> vector_intersect(std::vector<type> v1, std::vector<type> v2)
    {
        std::sort(v1.begin(), v1.end());
        std::sort(v2.begin(), v2.end());

        std::vector<type> v_intersection;

        std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), std::back_inserter(v_intersection));

        return v_intersection;
    }

    std::vector<MRRGNode *> filterCanMapOp(std::vector<MRRGNode *> fus, OpGraphOp const *op)
    {
        std::vector<MRRGNode *> result;
        for (auto &f : fus)
        {
            if (f->canMapOp(op))
            {
                result.push_back(f);
            }
        }

        return result;
    }
    
    void updateOveruse(MRRG *mrrg)
    {
        overuse_num_ = 0;
        for (auto &node : mrrg->routing_nodes)
        {
            if (occupancy[node] > node->capacity)
            {
                overuse_num_ += occupancy[node]-  node->capacity;
            }
        }
    }


protected:
    int rand_seed;
    float pfactor;
    float pfactor_factor;
    float const_temp_factor;
    int swap_factor;
    float cold_accept_rate;
    float updateTempConst(float temp);

private:
    bool inner_place_and_route_loop(OpGraph *opgraph, MRRG *mrrg, float temp, float *accept_rate);
    MRRGNode *getCandidateFU(MRRG *mrrg, OpGraphOp *op);
    MRRGNode *getRandomUnoccupiedFU(MRRG *mrrg, OpGraphOp *op);
    OpGraphOp *getOpNodePtr(OpGraph *opgraph, MRRGNode *n);

    OpMapping ripUpOp(OpGraphOp *op, float *cost = NULL);
    void restoreOp(OpMapping oldmap);

    bool routeOp(OpGraphOp *op, MRRG *mrrg);
    bool routeVal(OpGraphVal *val);
    bool placeOp(OpGraphOp *op, MRRGNode *n);

    bool checkOveruse(MRRG *mrrg);

    // mapping/unmapping
    void mapMRRGNode(OpGraphNode *, MRRGNode *node);
    void unmapMRRGNode(OpGraphNode *, MRRGNode *node);
    std::vector<MRRGNode *> unmapAllMRRGNodes(OpGraphNode *);
    void mapAllMRRGNodes(OpGraphNode *, std::vector<MRRGNode *> nodes);
    MRRGNode *getMappedMRRGNode(OpGraphOp *op);

    // mapping and occupancy
    std::map<MRRGNode *, int> occupancy;
    std::map<OpGraphNode *, std::vector<MRRGNode *>> mapping;

    // Costing
    float getTotalOpCost(OpGraphOp *op);
    float getCost(MRRGNode *n);
    float getCost(MRRG *n);
    float getCost(OpGraphNode *n);
    float getCost(OpGraph *opgraph);
    bool compare_mrrg_node_cost(MRRGNode *a, MRRGNode *b);

    MRRG * this_mrr;
    std::set<mpos> mapped_pos;

    bool allow_overuse = true;

    int MAX_OVERUSE = 1000; // I set 1000 as it is a good value to get an initial mapping.

    int overuse_num_ = 0;
};

#endif
