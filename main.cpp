//
// Created by caronkey on 10/5/2021.
//
#include "include/Graph.h"
#include "include/matching.h"
#include <iostream>
#include <fstream>
#include <map>
#include <cstring>
#include <sstream>
#include "include/cxxopts.hpp"
#include "include/Vertex.h"
int VERBOSE = 0;
const float ZERO = 0.00000001;
const float M_WEIGHT = 1000000;
float DEPTH = 100;
bool BREAK_C = false;
int TYPE = 0;
bool SELF_L = false;
int MIN_L = 1000;
int MODEL = 0;
bool DEBUG = false;
std::string SUB_ONLY = "";

std::unordered_map<std::string, int> max_seg_nodes(std::string &graph_file) {
    std::unordered_map<std::string, int> out_results;
    std::unordered_map<std::string, int> in_results;
    std::string source, target, startTag, originalSource, originalTarget;
    char sDir, tDir;
    int copyNum;
    float coverage = 1;
    float weight;
    float score;
    int cGene;
    std::ifstream infile(graph_file);
    std::string line;
    std::set<std::string> visited_self_loop;
    std::set<std::string> prev_junc;
    while (getline(infile, line)) {
        std::istringstream iss(line);
        if (line.rfind("SEG", 0) == 0) {
            iss >> startTag >> originalSource >> coverage >> copyNum >> cGene >> score;
            out_results[originalSource] = 0;
            in_results[originalSource] = 0;
        } else {
            iss>>startTag>>originalSource>>sDir>>originalTarget>>tDir>>weight;
            auto junc_str = originalSource+sDir+originalTarget+tDir;
            if (prev_junc.count(junc_str) != 0){
                continue;
            }
            prev_junc.insert(originalSource+sDir+originalTarget+tDir);
            if (sDir == '+') {
                if(out_results.find(originalSource) != out_results.end()) {
                    out_results[originalSource]++;
                }
                sDir = '-';
            } else {
                if(in_results.find(originalSource) != in_results.end()) {
                    in_results[originalSource]++;
                }
                sDir = '+';
            }
            if (tDir == '+') {
                if(in_results.find(originalTarget) != in_results.end()) {
                    in_results[originalTarget]++;
                }
                tDir = '-';
            } else {
                if(out_results.find(originalTarget) != out_results.end()) {
                    out_results[originalTarget]++;
                }
                tDir = '+';
            }
            prev_junc.insert(originalTarget+tDir+originalSource+sDir);
        }
    }
    for (auto& item : out_results) {
        auto k = item.first;
        auto v = item.second;
        if (k == "k141_19604") {
            auto tmp = 9;
        }
        if (in_results[k] > v) {
            item.second = in_results[k];
        }
    }
    return out_results;
}

bool check_visited_path(std::vector<std::vector<std::string>>& visited_vec, std::vector<std::string>& path){
    for (auto& p: visited_vec){
        if (p == path){
            return true;
        }
    }
    return false;
}

void checkMatrixConjugate(seqGraph::SparseMatrix& matrix, int n) {

//    for (auto i : matrix.IA) {
//        for (auto j : matrix.JA) {
//            auto t1 = matrix.getIJ(i,j);
//            auto t2 = matrix.getIJ(conjugateIdx(j),conjugateIdx(i));
//            if (std::abs(t1 - t2) > ZERO) {
//                std::cout<<"Error, not conjugated\n"<<i<<j<<std::endl;
//                exit(0);
//            }
//        }
//    }
    for(int i = 1; i < n+1; i++) {
        std::cout<<i<<" ";
        for(int j = 1; j < n+1; j++) {
            if(i == 3 && j == 15)
                auto oo = 33;
            auto t1 = matrix.getIJ(i,j);
            auto t2 = matrix.getIJ(seqGraph::conjugateIdx(j), seqGraph::conjugateIdx(i));
            if (std::abs(t1 - t2) > ZERO) {
                std::cout<<"Error "<<i<<" "<<j<<std::endl;
                exit(0);
            }
//            std::cout<<matrix[i][j]<<" ";
        }
        std::cout<<"\n"<<std::endl;
    }
}
void tokenize(const std::string &str, std::vector<std::string> &tokens, const std::string &delimiters)
{
    std::string::size_type pos, lastPos = 0, length = str.length();

    using value_type = typename std::vector<std::string>::value_type;
    using size_type  = typename std::vector<std::string>::size_type;

    while (lastPos < length + 1)
    {
        pos = str.find_first_of(delimiters, lastPos);
        if (pos == std::string::npos)
        {
            pos = length;
        }

        if (pos != lastPos)
            tokens.emplace_back(str.data() + lastPos, (size_type) pos - lastPos);

        lastPos = pos + 1;
    }
}
//recall paths from iteration result paths
//void pathsRecall(std::vector<std::map<int, std::vector<int>*>*> & recallPaths) {
//    for (auto iterPaths: *recallPaths.back()) {
//        for (auto item : *iterPaths.second) {
//
//        }
//    }
//}
//void print_help() {
//    std::cout<<"matching is a conjugated assembly algorithm that constructs a conjugate\n"
//               "map according to Barcode connection information between Contig,\n"
//               "and finds the most likely connection relationship between Contig\n";
//    std::cout<<"usage:\n\tmatching input.graph result.txt result.c.txt 0 1"<<std::endl;
//    std::cout<<"parameters:\n";
//    std::cout<<"\tinput.graph: graph generated by bDistance\n";
//    std::cout<<"\tresult.txt: linear output result\n";
//    std::cout<<"\tresult.c.txt: cycle output result\n";
//    std::cout<<"\t0: iteration matching times, 0 means only matching once, no more iteration\n";
//    std::cout<<"\t1: if output verbose information\n";
//}

void parse_tgs(std::string& f_name, seqGraph::Graph* g){
    std::ifstream tgs(f_name);
    std::string line;
    seqGraph::Vertex* vertex1 = nullptr;
    seqGraph::Vertex* vertex2 = nullptr;
    while( std::getline(tgs,line) )
    {
        std::stringstream ss(line);

        std::string v1, v2;
        char dir1, dir2;
        std::string course;
        int i = 1;
        while( std::getline(ss,course,' ') )
        {
            course.erase(std::remove(course.begin(), course.end(), '\r'), course.end());
            course.erase(std::remove(course.begin(), course.end(), '\n'), course.end());
            if (i == 1) {
                dir1 = course[course.size()-1];
                course.pop_back();
                v1 = course;
                vertex1 = g->addVertex(v1,"xx",1,2,1,1,1,0);
            } else {
                dir2 = course[course.size()-1];
                course.pop_back();
                v2 = course;
//                add to grap
                vertex2 = g->addVertex(v2,"xx",1,2,1,1,1,0);
                g->addJunction(vertex1, vertex2, dir1, dir2, M_WEIGHT, 1 , 1);
                vertex1 = vertex2;
                vertex2 = nullptr;
                dir1 = dir2;
//              make
            }
            i++;
        }
        std::cout<<"parse tgs down"<<"\n";
    }
    tgs.close();
}

int update_by_gene_junc_remote(seqGraph::Vertex& v1, seqGraph::Vertex& v2, char v1d, char v2d, seqGraph::Graph* g) {
    auto juncs_s = g->doesPathExists(v1,v2,v1d,v2d);
    bool is_updated = false;
    for (const auto& juncs : juncs_s) {
        float juncs_count = juncs.size();
        if (juncs_count == 0) return 0;
        float path_weight = DEPTH/juncs_count;
//    if (juncs_count > 5) return 0;
        for (auto junc: juncs) {
            junc->getWeight()->setCopyNum(junc->getWeight()->getCopyNum() + path_weight);
        }
        is_updated = true;
    }
    if(is_updated) return 1;
    return 0;
}

int parse_gene_junc(std::string& f_name, seqGraph::Graph* g, float theta) {
    std::string source, target, startTag, originalSource, originalTarget, prev_source;
    char sDir, tDir;
    int copyNum;
    float coverage = 1;
    int distance;
    float score;
    int cGene;
    std::ifstream infile(f_name);
//    this used for path backtrack
    std::string line;
    std::set<std::string> visited_self_loop;
    float distance_idx = 0;
    int is_updated = 0; //if A,B is updated we will skip all the next pair start from A
    while (getline(infile, line)) {
        std::istringstream iss(line);
        iss >> startTag >> originalSource >> sDir >> originalTarget >> tDir >> distance;
        if (originalSource == "H:960") {
            auto tmppp = 44;
        }
        source = originalSource;
        target = originalTarget;
        if (source == "EDGE_54893_length_473_cov_4.568182")
            auto mk = 33;
        auto v1t = g->getVertexById(source + "_0");
        auto v2t = g->getVertexById(target + "_0");
        if(v1t == nullptr or v2t == nullptr) continue;
//            v1 = v1t == nullptr ? g->addVertex(source,"xx",1,2,1,1,2) : v1t;
//            v2 = v2t == nullptr ? g->addVertex(target,"xx",1,2,1,1,2) : v2t;
        float c1 = v1t->getWeight()->getCopyNum();
        float c2 = v2t->getWeight()->getCopyNum();
        if (c1 == 0 || c2 == 0) {
            continue;
        }
        if (source == prev_source) {
            distance_idx++;
        } else {
            distance_idx = 1;
        }
//            divide with min copy and only assign
        auto min_copy = std::min(c1, c2);
//            auto divide_copy = std::min(c1,c2);
        auto multi_copy = c1 * c2;
        float gene_weight = DEPTH * (1/distance_idx) * theta;
        gene_weight = gene_weight / min_copy;
        source = originalSource;
        target = originalTarget;
//                if (MODEL == 2) {
// if the closest jun have gene
        if(source == "NODE_45_length_1665_cov_23.500000" && target == "NODE_46_length_1652_cov_28.206478_439_1652") {
            int tmp = 33;
        }
//        update_by_gene_junc_remote(*v1t,*v2t,sDir,tDir, g);
        if(!is_updated) {
            is_updated = update_by_gene_junc_remote(*v1t,*v2t,sDir,tDir, g);
        } else {
            if (prev_source != originalSource) {
                is_updated = update_by_gene_junc_remote(*v1t,*v2t,sDir,tDir, g);
            }
        }
        for (int i = 0; i < c1; i++) {
            v1t = g->getVertexById(source.append("_").append(std::to_string(i)));
//                g->addJunction(v1t, v2t, sDir, tDir, avg_weight, 1 , 1);
            for (int j = 0; j < c2; j++) {
                v2t = g->getVertexById(target.append("_").append(std::to_string(j)));
                target = originalTarget;
                auto jun = g->doesJunctionExist(*v1t, *v2t, sDir, tDir);
                if (jun == nullptr) continue;
                jun->getWeight()->setCopyNum(gene_weight+jun->getWeight()->getCopyNum());
            }
            source = originalSource;
        }
        prev_source = originalSource;
        // if a b c d, if a-d have a gene link,
    }
//                }
//            g->addJunction(v1t, v2t, sDir, tDir, weight, 1 , 1);
    infile.close();
}

//TODO add args support
// 1 graph 2 result 3 result.c 4 iteration 5 verbose 7 tgs
int main(int argc, char *argv[]) {

//    parse options
    cxxopts::Options options("matching", "Matching graph with conjugate bi-matching algorithm");

    options.add_options()
            ("g,graph", "Graph file", cxxopts::value<std::string>())
            ("r,result","Result", cxxopts::value<std::string>())
            ("c,result_c", "Cycle result", cxxopts::value<std::string>())
            ("i,iteration", "Iteration times", cxxopts::value<int>()->default_value("0"))
            ("v,verbose", "Verbose", cxxopts::value<int>()->default_value("0"))
            ("l,local_order", "Local order information", cxxopts::value<std::string>())
            ("gene_junc", "Gene junc", cxxopts::value<std::string>())
            ("m,result_m","Matching result", cxxopts::value<std::string>())
            ("b,break_c", "Whether break and merge cycle into line paths", cxxopts::value<bool>()->default_value("false"))
            ("s,self_l", "Cycle result", cxxopts::value<bool>()->default_value("false"))
            ("min_l", "Min length to print", cxxopts::value<int>()->default_value("-1"))
            ("sub_only", "Only get all sub graph",cxxopts::value<std::string>()->default_value(""))
            ("model", "0, meta; 1: single strain; ",cxxopts::value<int>()->default_value("0"))
            ("ignore_copy", "Ignore copy number", cxxopts::value<bool>()->default_value("false"))
            ("d,debug","debug user", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print usage");
    auto result = options.parse(argc,argv);
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    std::string graphF = result["graph"].as<std::string>();
    std::string resultF = result["result"].as<std::string>();
    std::string resultCF = result["result_c"].as<std::string>();
    int iterRounds = result["iteration"].as<int>();
    VERBOSE = result["verbose"].as<int>();
    MODEL = result["model"].as<int>();


    std::ifstream infile(graphF);
    std::ofstream resultFile(resultF);
    std::ofstream cyclePathsFile(resultCF);
    std::ofstream matchResultFile;
    int iterRoundsBK = iterRounds;
    if (result.count("result_m")) {
        std::string resultM = result["result_m"].as<std::string>();
        matchResultFile.open(resultM);
    }
    if (result.count("break_c")) {
        BREAK_C = true;
    }
    if (result.count("debug")) {
        DEBUG = true;
    }
    if (result.count("sub_only")) {
        SUB_ONLY = result["sub_only"].as<std::string>();
    }
    bool ignore_copy = false;
    if (result.count("ignore_copy")) {
        ignore_copy = result["ignore_copy"].as<bool>();
    }
    if (result.count("self_l")) {
        SELF_L = true;
    }
    auto* g = new seqGraph::Graph;
    if (result.count("local_order")) {
        std::string localOrder = result["local_order"].as<std::string>();
        parse_tgs(localOrder, g);
    }

    std::string source, target, startTag, originalSource, originalTarget;
    char sDir, tDir;
    int copyNum;
    float coverage = 1;
    float weight;
    float score;
    int cGene;

//
    auto seg_max_node = max_seg_nodes(graphF);
//    this used for path backtrack
    std::string line;
    std::set<std::string> visited_self_loop;
    while (getline(infile, line)) {
        std::istringstream iss(line);
        if (line.rfind("SEG", 0) == 0) {
//            print result match
            if (result.count("result_m")) {
                matchResultFile<<line<<std::endl;
            }
//            TODO, take coverage into consideration
            iss>>startTag>>originalSource>>coverage>>copyNum>>cGene>>score;
//            iss>>startTag>>originalSource>>copyNum;

//            TODO consider optimize the copied vertices.
            source = originalSource;
            source = source.append("_0");
            if(originalSource == "k141_19604") {
                int mm = 9;
            }
            auto copyNumNew = std::min(copyNum, seg_max_node[originalSource]);
            auto v = g->addVertex(source, "xx", 1, 2, coverage, 1, copyNumNew, 0, copyNum);
//            score = 1;
            v->setGeneAndScore(cGene, score);
            if (!ignore_copy) {
                for (int i = 1; i < copyNumNew; i++) {
                    source.pop_back();
                    if (i>=11)
                        source.pop_back();
                    if (i>=101)
                        source.pop_back();
                    source.append(std::to_string(i));
                    v = g->addVertex(source, "xx", 1, 2, coverage, 1, copyNumNew,i, copyNum);
                    v->setGeneAndScore(cGene, score);
//                source = originalSource;
                }
            }
        } else {
            iss>>startTag>>originalSource>>sDir>>originalTarget>>tDir>>weight;
            if (originalSource == "H:960") {
                auto tmppp = 44;
            }
            if(SELF_L && originalSource == originalTarget) {
                if (visited_self_loop.find(originalSource) == visited_self_loop.end()){
//                    std::vector<std::string> tokens;
//                    tokenize(originalSource, tokens, "_");
//                    std::cout<<originalSource<<std::endl;
//                    int length = std::stoi(tokens[3]);
//                    if (length>1000){
//                        std::cout<<"pass for: "<<originalSource<<std::endl;
//                        cyclePathsFile<< "self loop: "<<std::endl;
//                        cyclePathsFile<<originalSource+"+"<<std::endl;
//                    }else{
//                        std::cout<<"not pass for: "<<originalSource<<std::endl;
//                    }
                    visited_self_loop.insert(originalSource);
                }
                continue;
            }
//            seqGraph::Vertex* v1;
//            seqGraph::Vertex* v2;
            source = originalSource;
            target = originalTarget;
            if (source == "k141_19604")
                auto mk = 33;
            auto v1t = g->getVertexById(source+"_0");
            auto v2t = g->getVertexById(target+"_0");
//            v1 = v1t == nullptr ? g->addVertex(source,"xx",1,2,1,1,2) : v1t;
//            v2 = v2t == nullptr ? g->addVertex(target,"xx",1,2,1,1,2) : v2t;
//here, change copy to max seg copy
//            float c1 = v1t->getWeight()->getCopyNum();
//            float c2 = v2t->getWeight()->getCopyNum();
            float c1 = std::min(v1t->getWeight()->getCopyNum(), float(seg_max_node[originalSource]));
            float c2 = std::min(v2t->getWeight()->getCopyNum(), float(seg_max_node[originalTarget]));
            if (c1 == 0 || c2 == 0) {
                continue;
            }
//            divide with min copy and only assign
            auto min_copy = std::min(c1,c2);
//            auto divide_copy = std::min(c1,c2);
            auto multi_copy = c1 * c2;
            auto avg_weight = weight / min_copy;
            g->addJunction(v1t, v2t, sDir, tDir, avg_weight, 1 , 1);
            if (!ignore_copy) {
//                if (MODEL == 2) {
                    for (int i = 0; i < c1; i++) {
                        v1t = g->getVertexById(source.append("_").append(std::to_string(i)));
//                g->addJunction(v1t, v2t, sDir, tDir, avg_weight, 1 , 1);
                        for (int j = 0; j < c2; j++) {
                            v2t = g->getVertexById(target.append("_").append(std::to_string(j)));
                            g->addJunction(v1t, v2t, sDir, tDir, avg_weight, 1, 1);
                            target = originalTarget;
                        }
                        source = originalSource;
                    }
//                }
            }
//            g->addJunction(v1t, v2t, sDir, tDir, weight, 1 , 1);
        }
    }
    if (result.count("gene_junc")) {
        std::string gene_junc_file = result["gene_junc"].as<std::string>();
        parse_gene_junc(gene_junc_file, g, 1);
    }
//    matching for each connected graph
    int n = 0;
    g->parseConnectedComponents();
    if (SUB_ONLY != "") return 0;
    std::cout<<"total nodes"<<g->getVertices()->size()<<std::endl;
    int maxI = g->subGraphCount();
//    maxI = 3;
    std::cout<<"Isolated nodes"<<g->getIsolatedVs()->size()<<std::endl;
    for (auto item : *g->getIsolatedVs()) {
        resultFile<<item->getOriginId()<<"+"<<"\n";
    }
    std::vector<std::string> to_be_remove_sl;
    while (n < maxI) {
//        if(n == 118) continue;
//        if(n==6){
//            int m = 9;
//        }
        std::cout<<"process subgraph "<<n<<"\n";
        auto subGraph = g->getSubgraph(n);
//          auto subGraph = g;
//        subGraph->removeByGeneAndScore(cyclePathsFile);
        std::cout<<"sub graph nodes: "<<subGraph->getVertices()->size()<<std::endl;
        if(subGraph->getVertices()->size() == 1) {
            resultFile<<subGraph->getVertices()->front()->getOriginId()<<"+"<<"\n";
            n++;
            continue;
        }
//        if (subGraph->getJuncSize() == 1) {
//            auto sV = subGraph->getJunctions()->front()->getSource()->getOriginId();
//            auto tV = subGraph->getJunctions()->front()->getTarget()->getOriginId();
//            auto sD = subGraph->getJunctions()->front()->getSourceDir();
//            auto tD = subGraph->getJunctions()->front()->getTargetDir();
//            resultFile<<sV<<sD<<"\t"<<tV<<tD<<"\n";
//            n++;
//            continue;
//        }
        auto* m = new matching(subGraph);
//        m->checkConjugateMatrix();
//
        if (VERBOSE >= 2)
            checkMatrixConjugate(m->getMatrix(), m->getN());
//        if (g->getVertices()->size() >= 1000) {
//            m->main_steps();
//        } else {
//            m->hungarian();
//        }
        m->main_steps();
        if (VERBOSE >= 2) {
            std::cout<<"final matched relation\n";
            for(int i = 0; i < m->getN() + 1; i++) {
                std::cout<<m->getMatched()[i]<<"\t";
            }
        }
        std::cout<<"\nresolve path"<<std::endl;
//    TODO , make the path and cycles info into mathicng class
        if (result.count("result_m")) {
            m->writeMatchResult(matchResultFile);
        }
        auto paths = m->resolvePath(nullptr);

//        cyclePathsFile<<"sub\n";
//        cyclePathsFile<<"iter 0\n";
        std::set<std::string> seen;
        std::vector<std::vector<std::string>> visited_path;

        for (auto &item: *paths) {
            if (BREAK_C && m->isMergedPath(item.first)) continue;
            if (m->isCycle(item.first)) {
//                if (m->graph->getVertexByIdQ(m->idx2IdStr(item.second->front()))->sameVertex(*m->graph->getVertexByIdQ(m->idx2IdStr(item.second->back()))) ){
//                    item.second->pop_back();
//                }
                int idp = (&item-&*(paths->begin()));

                int total_length = 0;
                std::vector<int> path;

                bool not_print_tag;

                for(const auto& v: *item.second) {
//                    path.emplace_back(v)
                    int idx = (&v-&*(item.second->begin()));
                    // record IDX OF A-B-A_0-B_0
//                    if (m->idx2StrDir(v)=="EDGE_74036_length_196_cov_93.191489+"){
//                        std::cout<<'fk'<<std::endl;
//                        if (idx!=0 ){
//                            std::cout<<m->idx2StrDir(item.second->front())<<std::endl;
//                            std::cout<<m->idx2StrDir(v)<<std::endl;
//                            std::cout<<m->idx2StrDir((*item.second)[idx-1])<<std::endl;
//                            std::cout<<m->idx2StrDir(item.second->back())<<std::endl;
//                        }
//                    }
                    std::string true_name_front;
                    std::string true_name_current;
                    std::string true_name_pre;
                    std::string true_name_back;

                    if(idx!=0){
                        std::string idStrFront = m->idx2Str(item.second->front());
                        auto pos_front = idStrFront.find_last_of('_');
                        true_name_front = idStrFront.substr(0,pos_front);

                        std::string idStrCurrent = m->idx2Str(v);
                        auto pos_current = idStrCurrent.find_last_of('_');
                        true_name_current = idStrCurrent.substr(0,pos_current);

                        std::string idStrPre = m->idx2Str((*item.second)[idx-1]);
                        auto pos_pre = idStrPre.find_last_of('_');
                        true_name_pre = idStrPre.substr(0,pos_pre);

                        std::string idStrBack = m->idx2Str(item.second->back());
                        auto pos_back = idStrBack.find_last_of('_');
                        true_name_back = idStrBack.substr(0,pos_back);
                    }


                    if (idx!=0 && true_name_front==true_name_current && true_name_pre == true_name_back){
                        continue;
                        not_print_tag = true;
                    }
                    else{
                        not_print_tag = false;
                    }


//
//                    if (idx!=0 && m->idx2VertexInCurrentGraph(item.second->front())->sameVertex(*(m->idx2VertexInCurrentGraph(v))) && m->idx2VertexInCurrentGraph((*item.second)[idx-1])->sameVertex(*(m->idx2VertexInCurrentGraph(item.second->back())))){
//                        not_print_tag = true;
//                    }
//                    else{
//                        not_print_tag = false;
//                    }
                    if (MODEL == 0) {
                        std::vector<std::string> tokens;
                        tokenize(m->idx2StrDir(v).substr(0, m->idx2StrDir(v).length()-1), tokens, "_");
                        int length = std::stoi(tokens[3]);
                        total_length+=length;
                    }
                }
                if(MODEL ==0 and (total_length<MIN_L|| not_print_tag)){
                    continue;
                }
//                bool self_loop_flag = false;
//                std::string path_str;
                std::vector<std::string> tmp;
                for(auto &t: *item.second){
                    std::cout<<m->idx2Str(t)<<std::endl;
                    std::string idStr = m->idx2Str(t);
                    auto pos = idStr.find_last_of('_');
                    std::string true_name = idStr.substr(0,pos);
                    tmp.emplace_back(true_name);
                }
                for(const auto& v: *item.second) {
//                    path_str+=("_"+std::to_string(v));
                    std::cout<<m->idx2Str(v)<<std::endl;
                    std::string idStr = m->idx2Str(v);
                    auto pos = idStr.find_last_of('_');
                    std::string true_name = idStr.substr(0,pos);
                    if (visited_self_loop.find(true_name)!=visited_self_loop.end()){
                        to_be_remove_sl.emplace_back(true_name);
//                        self_loop_flag = true;
                    }
                }
//                if(self_loop_flag) continue;
                std::string res_str;
//                if (std::find(visited_path.begin(), visited_path.end(), item.second) == visited_path.end()){
//                    cyclePathsFile<<"iter "<<0<<",graph"<<n<<"\n";
                if (!check_visited_path(visited_path, tmp)){
                    for(const auto& v: *item.second) {
                        int idr = (&v-&*(item.second->begin()));
//                        if (std::find(not_print_idxes.begin(), not_print_idxes.end(), idr) == not_print_idxes.end()) continue;
//                        cyclePathsFile<<m->idx2StrDir(v)<<"\t";
                        res_str+=(m->idx2StrDir(v)+"\t");
                    }
//                    cyclePathsFile<<"\n";
                    if(res_str.length()>0){
                        cyclePathsFile<<"iter "<<0<<",graph"<<n<<"\n";
                        cyclePathsFile<<res_str<<"\n";
                    }

                    std::sort(tmp.begin(), tmp.end());
                    visited_path.emplace_back(tmp);


               }

            }
        }

//    recallPaths.push_back(paths);
        int iterN = 0;
        iterRounds = iterRoundsBK;
        int prevPathSize = 0;
        if (paths->size() != 1) {
            while (iterRounds !=0) {
                std::cout<<"Iteration "<<iterN + 1<<",nodes count"<<paths->size()<<std::endl;
//                TODO if v==1 or juncs ==1 continue
                m->reconstructMatrix(paths, subGraph);
                if (m->getGraph()->getJuncSize() == 0) {
                    std::cout << "No match needed" << std::endl;
                    break;
                }else{
                    //                checkMatrixConjugate(m->getMatrix(), m->getN());
                    std::cout<<"start hungarian"<<std::endl;
                    //                m->hungarian();
                    m->main_steps();
                }
                paths = m->resolvePath(paths);
                if (paths->size() == prevPathSize || paths->size() == 1) {std::cout<<paths->begin()->second->size()<<"break"<<std::endl; break;}
                for (auto item: *paths) {
                    if (BREAK_C && m->isMergedPath(item.first)) continue;
                    if (m->isCycle(item.first)) {
                        cyclePathsFile<<"iter "<<iterN<<",graph"<<n<<"\n";
                        for(const auto& v: *item.second) {
                            cyclePathsFile << m->idx2StrDir(v) << "\t";
                        }
                        cyclePathsFile<<"\n";
                    }
                }
                prevPathSize = paths->size();
                iterRounds--;
                iterN++;
            }
        }
        std::vector<std::string> all_paths;
        for(auto item : *paths) {
            if (BREAK_C && m->isMergedPath(item.first)) continue;
            int total_length = 0;
            std::string current_path;
//            if copy = 1 or copy > and nextjunc == 1, add, else break
            for (int i = 0; i < item.second->size(); i ++) {
                auto v = (*item.second)[i];
//                auto v_next = (*item.second)[i + 1];
                if (MODEL == 0) {
                    std::vector<std::string> tokens;
                    tokenize(m->idx2StrDir(v).substr(0, m->idx2StrDir(v).length()-1), tokens, "_");
                    int length = std::stoi(tokens[3]);
                    total_length+=length;
                }
                if (v == -1) {
                    continue;
                }
//                if (MIN_L != -1 && total_length < MIN_L) continue;
                int v_copy = m->idx2VertexInOriginalGraph(v)->getWeight()->getCopyNum();
                auto tmp = m->idx2StrDir(v);

                if (tmp == "NZ_CP049085.2:1-114526+") {
                    int ttt = 0;
                }
                int v_outdegree = m->outDegree(v);
                int v_indegree = m->inDegree(v);
                if (DEBUG && item.second->size() > 2) {
//                    if ((v_copy == 1 || (v_copy > 1 && (v_outdegree <= 1 && v_indegree <= 1)))) {
//                    if (((v_outdegree <= 1 && v_indegree <= 1))) {
//                        current_path.append(m->idx2StrDir(v)).append("\t");
//                    } else {
//                        current_path.append("\n");
//                        current_path.append(m->idx2StrDir(v)).append("\t");
//                    }
                    if (tmp == "NODE_51_length_1238_cov_40.515935_96_1238+") {
                        int ttt = 0;
                    }
                    bool isprinted = false;
                    if ((v_indegree > 1)) {
                        current_path.append("\n");
                        current_path.append(tmp).append("\t");
                        isprinted = true;
                    }
                    if (((v_outdegree > 1))) {
                        if(!isprinted) {
                            current_path.append(tmp).append("\t");
                        }
                        current_path.append("\n");
                    } else {
                        if(!isprinted)
                            current_path.append(m->idx2StrDir(v)).append("\t");
                    }

                } else {
                    current_path.append(m->idx2StrDir(v)).append("\t");
                }
            }
//            for(const auto& v: *item.second) {
//                if (MODEL == 0) {
//                    std::vector<std::string> tokens;
//                    tokenize(m->idx2StrDir(v).substr(0, m->idx2StrDir(v).length()-1), tokens, "_");
//                    int length = std::stoi(tokens[3]);
//                    total_length+=length;
//                }
//                if (v == -1) {
//                    resultFile<<'c';
//                    continue;
//                }
////                if (MIN_L != -1 && total_length < MIN_L) continue;
//                int v_indegree = m->inDegree(v);
//                int v_outdegree = m->outDegree(v);
//                if (DEBUG && v_indegree > 1) {
//                    current_path.append("\n");
//                }
//                current_path.append(m->idx2StrDir(v)).append("\t");
//                if (DEBUG && v_outdegree > 1) {
//                    current_path.append("\n");
//                }
//            }
            if (MODEL == 0) {
                if (total_length>500){
                    if (std::find(all_paths.begin(), all_paths.end(), current_path) != all_paths.end())
                        continue;
                    all_paths.push_back(current_path);
                    resultFile<<current_path<<std::endl;
                }else{
                    std::cout<<"not pass for: "<<current_path<<std::endl;
                }
            } else {
                if (std::find(all_paths.begin(), all_paths.end(), current_path) != all_paths.end())
                    continue;
                all_paths.push_back(current_path);
                resultFile<<current_path<<std::endl;
            }
        }
        free(m);
        free(subGraph);
        n++;
    }
    for(auto& self_loop: visited_self_loop){
        if (std::find(to_be_remove_sl.begin(), to_be_remove_sl.end(), self_loop)==to_be_remove_sl.end()){
            if (MODEL == 0) {
                std::vector<std::string> tokens;
                tokenize(self_loop, tokens, "_");
                std::cout<<self_loop<<std::endl;
                int length = std::stoi(tokens[3]);
                if (length>1000){
                    std::cout<<"pass for: "<<self_loop<<std::endl;
                    cyclePathsFile<< "self loop: "<<std::endl;
                    cyclePathsFile<<self_loop+"+"<<std::endl;
                }else{
                    std::cout<<"not pass for: "<<self_loop<<std::endl;
                }
            } else {
                cyclePathsFile<< "self loop: "<<std::endl;
                cyclePathsFile<<self_loop+"+"<<std::endl;
            }
        }
    }
    infile.close();
    resultFile.close();
    cyclePathsFile.close();
    matchResultFile.close();
}
