#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <queue>
#include <set>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "MeshAdjReorder.h"

bool hasNode(apf::Mesh* m, apf::MeshEntity* e)
{
    return m->getShape()->countNodesOn(m->getType(e)) > 0;
}

void getGoodStartingElement(apf::Mesh* m, apf::MeshEntity* & start, int dim) {
    apf::MeshIterator* it = m->begin(dim);
    /*below removed for cases where geometric model does not
    * have any vertices, because the below depends on that invariant
    * so we disable it until later when our code guarentees it*/
    #if 0
    // unsigned int min_order = 400;
    // while((e = m->iterate(it))){
    //     apf::ModelEntity* model_ent = m->toModel(e);
    //     assert(model_ent != NULL);
    //     int dim = m->getModelType(model_ent);
    //     if(dim != 0) {
    //         std::cout << "continue" << std::endl;
    //         continue;
    //     }
    //     apf::Adjacent adj;
    //     m->getAdjacent(e, 1, adj);
    //     //this keeps the first lowest order is sees
    //     std::cout << adj.getSize() << std::endl;
    //     if(adj.getSize() < min_order) {
    //         min_order =  adj.getSize();
    //         start= e; 
    //     }
    // }
    #endif
    /* return the first vertex*/
    start = m->iterate(it);
    m->end(it);
    // std::cout << start << std::endl;
}

struct sortVertsFunctor {
    apf::Mesh* m;
    bool operator()(apf::MeshEntity *lhs, apf::MeshEntity *rhs) {
    apf::Adjacent adj;
    this->m->getAdjacent(lhs, 1, adj);
    int lhs_count = adj.getSize();
    this->m->getAdjacent(rhs, 1, adj);
    int rhs_count = adj.getSize();
    return (lhs_count < rhs_count);
    }
};

void adjReorder(apf::Mesh* m, 
            apf::FieldShape* shape, 
            uint32_t n_components,
            apf::Numbering* & nodeNums,
            apf::Numbering* & faceNums) 
{
    apf::MeshEntity* e;
    apf::MeshIterator* it;
    signed long node_label = 0;
    int dim = m->getDimension();

    for(int ii = 0; ii <= dim; ++ii){
        if(shape->hasNodesIn(ii)){
            it = m->begin(ii);
            while((e = m->iterate(it))) {
                // std::cout << e << std::endl;
                if(m->isOwned(e)) {
                    node_label += shape->countNodesOn(m->getType(e)) * n_components;
                }
            }
            m->end(it);
        }
    }
    

    /* deduct any fixed entities from the locally-owned total */
    apf::MeshEntity* start_vert = NULL;
    getGoodStartingElement(m, start_vert, 0);
    // std::cout << start_vert << std::endl;
    std::queue<apf::MeshEntity*> q;
    std::set<apf::MeshEntity*> s;
    q.push(start_vert);
    s.insert(start_vert);
    std::set<apf::MeshEntity*>::iterator set_it;
    std::vector<apf::MeshEntity*> ent_list;

    /*switch to zero based indexing*/
    node_label -= 1;
    signed int face_label = m->count(2) - 1;
    /*
    * we stop at zero because we labeled every node
    * everything left in the queue will already have been labeled
    */
    while( q.size() > 0 && node_label >= 0) {
        apf::MeshEntity* current_entity;
        current_entity = q.front();
        q.pop();
        set_it = s.find(current_entity);
        if(set_it != s.end()) {
            s.erase(set_it);
        }
        /*assuming that there is only one dof per entity*/
        if(!apf::isNumbered(nodeNums, current_entity,0,0)) {
            for(uint32_t kk = 0; kk < n_components; ++kk){
                apf::number(nodeNums, current_entity, 0, (n_components - kk - 1), node_label);
                --node_label;
            }
        } else {
            continue;
        }

        /*
        * now we want to add any unlabled adjacent mesh entities to the queue
        * additions to queue will only be handled  by looking at adjacencies from
        * vertices
        */
        if(apf::getDimension(m, current_entity) == 0) {
            apf::MeshEntity* vert = current_entity;
            apf::Adjacent edge_adj;
            m->getAdjacent(vert, 1, edge_adj);
            int edge_adj_size = edge_adj.getSize();
            for(int ii = 0; ii < edge_adj_size; ++ii) {
                apf::MeshEntity* curr_edge = edge_adj[ii];
                //get the face adjacency for this edge and see if it is labeled
                apf::Adjacent face_adj;
                m->getAdjacent(edge_adj[ii], 2, face_adj);
                int face_adj_size = face_adj.getSize();
                for(int jj = 0; jj < face_adj_size; ++jj) {
                    if(!apf::isNumbered(faceNums, face_adj[jj], 0, 0)) {
                    apf::number(faceNums, face_adj[jj], 0 ,0, face_label);
                    --face_label;
                    }
                    //check if the face has a node to label, then queue it if exists
                    if(hasNode(m, face_adj[jj])) { // only for lagrange elements
                    //queue the face node
                    ent_list.push_back(face_adj[jj]);
                    //std::cout << "trying to queue a face node" << std::endl;
                    }
                }
                /*get the other end of the edge */
                apf::MeshEntity* other_vert = apf::getEdgeVertOppositeVert(m, curr_edge, vert);
                if(hasNode(m, curr_edge)) {
                    /*check that the other_vert is labeled or in queue and edge node not labeled*/
                    set_it = s.find(other_vert);
                    bool is_s = (set_it != s.end());
                    if( (apf::isNumbered(nodeNums, other_vert, 0, 0) || is_s ) &&
                        !apf::isNumbered(nodeNums, curr_edge, 0, 0) ) {
                        for(uint32_t kk = 0; kk < n_components; ++kk) {
                            apf::number(nodeNums, curr_edge, 0 ,  (n_components - kk - 1), node_label);
                            --node_label;
                        }
                    } else {
                        /*check that this edge is not labeled or in queue*/
                        set_it = s.find(curr_edge);
                        if(!apf::isNumbered(nodeNums, curr_edge, 0, 0) && set_it == s.end()) {
                            q.push(curr_edge);
                            s.insert(curr_edge);
                        }
                        /*add the other vertex to the list*/
                        ent_list.push_back(other_vert);
                    }
                } else {
                    if(!apf::isNumbered(nodeNums, other_vert, 0, 0)) {
                        ent_list.push_back(other_vert);
                    }
                }
            }/* finished loop over all edges coming into vertex */
            /* add all of the elements in the list to the queue */
            if(ent_list.size() > 0) {
                /*use the functor to create a closure over the mesh for ranking the vector */
                struct sortVertsFunctor functor;
                functor.m = m;
                std::sort(ent_list.begin(), ent_list.end(), functor);
                /*add each element of the vector to the queue*/
                for(std::vector<apf::MeshEntity*>::iterator it = ent_list.begin(); it != ent_list.end(); ++it) {
                    /*check that entity is not already in queue*/
                    set_it = s.find(*it);
                    if(set_it == s.end()) {
                        q.push(*it);
                        s.insert(*it);
                        } else {
                    // std::cout << "already in queue" << std::endl;
                    }
                }
                ent_list.clear();
            }
        } else {
            // std::cout << "not a vertex" << std::endl;
        } /* end of vertex checking*/
        if(q.size() > 500 ) {
            // std::cout << "bad termination" << std::endl;
            break;
        }
    } /*end of while loop */
    /*make sure that we did not leave an unlabeled node in the queue */
    while(!q.empty()) {
        for(uint32_t kk = 0; kk < n_components; ++kk){
            assert(apf::isNumbered(nodeNums, q.front(), 0, kk) == 1);
        }
        q.pop();
    }
}
