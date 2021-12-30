//
// Created by caronkey on 10/5/2021.
//
#include "include/Graph.h"
#include "include/matching.h"
#include <iostream>
#include <fstream>
#include <map>
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
            tokens.push_back(value_type(str.data() + lastPos, (size_type) pos - lastPos));

        lastPos = pos + 1;
    }
}

int main(int argc, char *argv[]) {
    auto md = argv[1];
    std::ifstream infile(argv[1]);
    std::string source, target;
    char sDir, tDir;
    double weight;
    seqGraph::Graph* g = new seqGraph::Graph;

    auto* idToSeg = new std::map<int, std::string>();
    auto* segToId = new std::map<std::string, int>();

    int idx = 0;
    while (infile>>source>>sDir>>target>>tDir>>weight) {
        seqGraph::Vertex* v1;
        seqGraph::Vertex* v2;
        auto v1t = g->getVertexById(source);
        auto v2t = g->getVertexById(target);
        v1 = v1t == nullptr ? g->addVertex(source,"xx",1,2,1,1,2) : v1t;
        v2 = v2t == nullptr ? g->addVertex(target,"xx",1,2,1,1,2) : v2t;
        g->addJunction(v1, v2, sDir, tDir, weight, 1 , 1);
    }
//    auto v1 = g->addVertex(1,"xx",1,2,1,1,2);
//    auto v2 =g->addVertex(2,"xx",1,2,1,1,2);
//    auto v3 =g->addVertex(3,"xx",1,2,1,1,2);
//    auto v4 =g->addVertex(4,"xx",1,2,1,1,2);
//    auto v5 =g->addVertex(4,"xx",1,2,1,1,2);
//
//    g->addJunction(v1,v2,'+','+',2,1,1);
//    g->addJunction(v2,v3,'+','-',2,1,1);
//    g->addJunction(v3,v4,'-','+',2,1,1);
//    g->addJunction(v4,v5,'+','-',2,1,1);

    matching* m = new matching(g);
    for(int i = 0; i < m->getN() + 1; i++){
        for(int j = 0; j < m->getN() + 1; j++) {
            auto tmp = m->getMatrix()[i][j];
            std::cout<<tmp<<"\t";
        }
        std::cout<<"\n";
    }
    std::cout<<std::endl;
    for(int i = 0; i < m->getN() + 1; i++) {
        std::cout<<m->getMatched()[i]<<"\t";
    }
    std::cout<<std::endl;

    m->main_steps();
    std::cout<<"out\n";
    for(int i = 0; i < m->getN() + 1; i++) {
        std::cout<<m->getMatched()[i]<<"\t";
    }
    std::cout<<"resolve\n";
    m->resolvePath();
}