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

#ifndef ___LISAMapper_H__
#define ___LISAMapper_H__

#include <string>
#include <map>
#include <vector>
#include <memory>
#include <algorithm>
#include <chrono>
#include <random>
#include <cmath>
#include <ctime>

#include <CGRA/CGRA.h>
#include <CGRA/OpGraph.h>
#include <CGRA/Mapper.h>
#include <CGRA/Mapping.h>
#include <CGRA/lisa/gnn.h>
#include <CGRA/lisa/LISAController.h>

using namespace std;

struct op_edge{
    OpGraphOp * src;
    OpGraphOp * des;
    int length;

   
};

struct op_edge_comp {
    bool operator()(op_edge const & a, op_edge const & b) {
        // your code here
        if(a.length != b.length){
            return a.length >  b.length;
        }else{
            // just give an result. 
            std::string as = a.src->name + a.des->name;
            std::string bs = b.src->name + b.des->name;
            return as.compare(bs);
        }
    }
};

class Generator {
    std::default_random_engine generator;
    std::normal_distribution<double> distribution;
    double min;
    double max;
public:
    Generator(double mean, double stddev, double min, double max):
        distribution(mean, stddev), min(min), max(max)
    {}

    double operator ()() {
        while (true) {
            double number = this->distribution(generator);
            if (number >= this->min && number <= this->max)
                return number;
        }
    }
};


class LISAMapper : public Mapper
{
public:
    // Custom setup with variable args
    LISAMapper(std::shared_ptr<CGRA> cgra, int timelimit, const std::map<std::string, std::string> &args);
    // Default setup
    LISAMapper(std::shared_ptr<CGRA> cgra);
    // Custom setup
    LISAMapper(std::shared_ptr<CGRA> cgra, int timelimit, int rand_seed, float initial_penalty, float penalty_factor, float const_temp_factor, int swap_factor, float cold_accept_rate);
    ~LISAMapper() = default;

    Mapping mapOpGraph(std::shared_ptr<OpGraph> opgraph, int II, std::string arch_model_name) override;


    float getFinalCost();

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

    void initLisa(std::vector<std::pair<int, int>> lisa_edges, std::set<std::pair<int, int>> lisa_backedges, std::map<OpGraphOp*, int>  node_to_id,  std::map<int, OpGraphOp*>  id_to_node, std::shared_ptr<std::map<int, node_label>> dfg_label){
        dfg_label_ = dfg_label;
        lisa_backedges_ = lisa_backedges;
        lisa_edges_ = lisa_edges;
        node_to_id_ = node_to_id;
        id_to_node_ = id_to_node;

        for(auto node: node_to_id_){
            node_children_[node.second] =  std::set<int>();
        }
        for(auto edge: lisa_edges){
            node_children_[edge.first].insert(edge.second);
        }
       
       
    }

    template <typename type>
    std::vector<type> vector_intersect(std::vector<type> v1, std::vector<type> v2);
    MRRGNode *  getRandomFU(MRRG *mrrg, OpGraphOp *op);
    // MRRGNode *  getRandomFUWithII(MRRG *mrrg, OpGraphOp *op);
    
    float updateTemperature(float t, float acceptance_rate);
    void topological_sort(OpGraph *g);
    std::vector<std::pair<MRRGNode *, int>> candidate_fu_intersect(std::vector<std::pair<MRRGNode *, int>> v1, std::vector<std::pair<MRRGNode *, int>> v2);
    std::vector<MRRGNode *> filterCanMapOp(std::vector<MRRGNode *> fus, OpGraphOp const *op);

    std::map<int, pos3d> dumpMapping(std::shared_ptr<OpGraph> opgraph, int &max_latency);
    std::string printMapping(){
        int test_latency =0 ;
        auto dumpedmapping = dumpMapping(opgraph_, test_latency);
        // std::cout<<"dumpedmapping "<<dumpedmapping.size()<<"\n";
        std::stringstream result;
        for(auto node: opgraph_->op_nodes){
            if(mapping.find(node)!=mapping.end() && mapping[node].size() > 0){
                result<<node->name<<" "<<mapping[node].front()->getFullName()<<"\n";
            }
            
        }
         result<< lisa_ctrl->mappingToStr(dumpedmapping, test_latency);
         return result.str();
    }

    bool optimizeMapping(OpGraph *opgraph, MRRG *mrrg);
    MRRGNode *  getLISAFU(std::shared_ptr<OpGraph> opgraph, MRRG *mrrg, OpGraphOp *op, int accepted = 1 , int total_tried =1, int num_swap = 1);
    std::vector<MRRGNode*> getDesiredFu(MRRG *mrrg, int start_II, int end_II, OpGraphOp *op);
    std::pair<int,int> getIntervalByScheduleOrder(std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, OpGraphOp *op, int scheduler_order );

    std::map<MRRGNode *, int> getCostByComm(MRRG *mrrg, std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op );
    std::map<MRRGNode *, int> getCostByAssociation(std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op, int start_II );
    std::map<MRRGNode *, int> getCostForSameLevelNode(std::shared_ptr<OpGraph> opgraph, std::map<int, pos3d> & dumped_mapping, std::vector<MRRGNode *> candidates, OpGraphOp *op );

     MRRGNode *  getCloseRandomFU(MRRG *mrrg, OpGraphOp *op, MRRGNode * old_fu){
        return  getCloseRandomFU(mrrg, op,  old_fu,  cgra_x_,  cgra_y_ );
     }
    MRRGNode *  getCloseRandomFU(MRRG *mrrg, OpGraphOp *op, MRRGNode * old_fu, int max_physical_dis,  int max_temp_dis );
     MRRGNode* getRoutingNode(MRRG *mrrg, int x, int y, int t);

    bool enable_eval_routing_priority(){ lisa_eval_routing_priority = true;}


    bool routeOp_withfailed_node(OpGraphOp *op, MRRG *mrrg, std::set<OpGraphOp *> &failed);

    void setLISAController( std::shared_ptr<LISAController> ctrl,int cgra_x,int cgra_y){ 
        lisa_ctrl = ctrl;
        cgra_x_ = cgra_x;
        cgra_y_ = cgra_y;
         node_asap_ = lisa_ctrl->getNodeASAP();
        }
    void enableTraining(){ is_training = true;}
    void disableTraining(){ is_training = false;}

    int MOD(int cycle){
        return (cycle + 100 * II_)%II_;
    }

     bool place_and_routing_parallel =  false;

    
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
    bool routeOpInput(OpGraphOp *op, MRRG *mrrg);
    bool routeVal(OpGraphVal *val,  OpGraphOp *output_op = NULL) ;
    bool placeOp(OpGraphOp *op, MRRGNode *n);

    bool checkOveruse(MRRG *mrrg);

    // mapping/unmapping
    void mapMRRGNode(OpGraphNode *, MRRGNode *node);
    void unmapMRRGNode(OpGraphNode *, MRRGNode *node);
    std::vector<MRRGNode *> unmapAllMRRGNodes(OpGraphNode *);
    void mapAllMRRGNodes(OpGraphNode *, std::vector<MRRGNode *> nodes);
    MRRGNode *getMappedMRRGNode(OpGraphOp *op);

    // mapping and occupancy
    // std::map<MRRGNode *, std::vector<OpGraphNode *>> occupancy_detail;
    std::map<MRRGNode *, int> occupancy;
    std::map<OpGraphNode *, std::vector<MRRGNode *>> mapping;

    // Costing
    float getTotalOpCost(OpGraphOp *op);
    float getCost(MRRGNode *n);
    float getCost(MRRG *n);
    float getCost(OpGraphNode *n);
    float getCost(OpGraph *opgraph);
    bool compare_mrrg_node_cost(MRRGNode *a, MRRGNode *b);



    //for lisa
    std::vector<std::pair<int, int>> lisa_edges_;//useless
    std::set<std::pair<int, int>> lisa_backedges_;//useless
    std::map<OpGraphOp*, int>  node_to_id_;
    std::map<int, OpGraphOp*>  id_to_node_;
    std::map<int, int> node_asap_;

    std::map<OpGraphOp*, int> overuse_counter;


    std::map<int, std::set<int>> node_children_;

    std::shared_ptr<std::map<int, node_label>> dfg_label_;

    std::set<op_edge, op_edge_comp> op_edges_;

    bool allow_overuse = false;

    int MAX_OVERUSE = 1000;

    int overuse_num_ = 0;

    int II_ = 0;

    MRRG *mrrg_ = NULL;

    const int MAX_SECOND_IN_ROUTING = 2;

    std::set<MRRGNode *> overuse_fu;

    std::set<mpos> mapped_pos;

    std::shared_ptr<LISAController> lisa_ctrl;

    int cgra_x_ = 0;
    int cgra_y_ = 0;

    std::shared_ptr<OpGraph> opgraph_;

    bool is_training = false;

    bool finish_init = false;


    bool lisa_eval_routing_priority = false;

   

    std::chrono::steady_clock::time_point mapper_start_time;

     std::set<OpGraphOpCode> lat_opcode ={OpGraphOpCode::OPGRAPH_OP_NOP ,
    OpGraphOpCode::OPGRAPH_OP_SEXT,
    OpGraphOpCode::OPGRAPH_OP_ZEXT,
    OpGraphOpCode::OPGRAPH_OP_TRUNC,
    OpGraphOpCode::OPGRAPH_OP_PHI,
    OpGraphOpCode::OPGRAPH_OP_ADD,
    OpGraphOpCode::OPGRAPH_OP_SUB,
    OpGraphOpCode::OPGRAPH_OP_MUL,
    OpGraphOpCode::OPGRAPH_OP_DIV,
    OpGraphOpCode::OPGRAPH_OP_AND,
    OpGraphOpCode::OPGRAPH_OP_OR,
    OpGraphOpCode::OPGRAPH_OP_XOR,
    OpGraphOpCode::OPGRAPH_OP_SHL,
    OpGraphOpCode::OPGRAPH_OP_SHRA,
    OpGraphOpCode::OPGRAPH_OP_SHRL,
    OpGraphOpCode::OPGRAPH_OP_LOAD,
    OpGraphOpCode::OPGRAPH_OP_STORE,
    OpGraphOpCode::OPGRAPH_OP_GEP,
    OpGraphOpCode::OPGRAPH_OP_ICMP,
    } ;


};

#endif
