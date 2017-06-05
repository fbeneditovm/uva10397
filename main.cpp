#include <stdio.h>
#include <vector>
#include <string>
#include "boost/tuple/tuple.hpp"

#include <limits.h>

#include <queue>
#include <algorithm>
#include <iostream>
#include <math.h>

using namespace std;
using boost::tuple;

#define MAXV 201

#define __DEBUG__

typedef struct {
    int V, E;
    int dg[MAXV];
    int adj[MAXV][MAXV], weight[MAXV][MAXV];
} Graph;

void initGraph(Graph *g, int v) {
    g->V = v;
    g->E = 0;
    for(int i=0; i<MAXV; i++)
        g->dg[i] = 0;
}

void addEdge(Graph *g, int u, int v, int w=1) {
    g->E++;
    g->adj[u][ g->dg[u] ] = v;
    g->weight[u][ g->dg[u] ] = w;
    g->dg[u]++;
}

void bfs(Graph *g, int node) {
    queue<int> q;
    int visited[MAXV];
    for(int i=0; i<MAXV; i++) {
        visited[i] = 0;
    }
    q.push(node);
    visited[node] = 1;

    while(!q.empty()) {
        int w = q.front();
        q.pop();
#ifdef __DEBUG__
        printf("BFS(%d)\n", w);
#endif
        for(int i=0; i<g->dg[w]; i++) {
            int v = g->adj[w][i];
            if(!visited[v]) {
                visited[v] = 1;
                q.push(v);
            }
        }
    }
}

void dfs_visit(Graph *g, int visited[MAXV], int node) {
#ifdef __DEBUG__
    printf("DFS(%d)\n", node);
#endif

    visited[node] = 1;
    for(int i=0; i<(g->dg[node]); i++) {
        int v = g->adj[node][i];
        if(!visited[v])
            dfs_visit(g, visited, v);
    }
}

void dfs(Graph *g, int node) {
    int visited[MAXV];
    for(int i=0; i<MAXV; i++) {
        visited[i] = 0;
    }
    dfs_visit(g, visited, node);
}

void dijkstra(Graph *g, int node, int cost[MAXV]) {
    int visited[MAXV];
    for(int i=0; i<MAXV; i++) {
        cost[i] = INT_MAX;
        visited[i] = 0;
    }
    cost[node] = 0;

    priority_queue< pair<int,int>,
            vector< pair<int,int> >,
            greater< pair<int,int> > >
            q;
    q.push( pair<int,int>(0, node) );

    while(!q.empty()) {
        //escolhe node de trabalho w
        /*int w = -1;
        for(int i=1; i<=g->V; i++) {
            if(!visited[i]) {
                if(w==-1)
                    w = i;
                else if(cost[i]<cost[w])
                    w = i;
            }
        }
        if(w==-1) break;*/

        pair<int,int> p = q.top();
        q.pop();
        int w = p.second;
        if(visited[w]) continue;

        visited[w] = 1;
        for(int i=0; i<g->dg[w]; i++) {
            int v = g->adj[w][i];
            int weightEdge = g->weight[w][i];
            if(cost[w]+weightEdge<cost[v]) {
                cost[v] = cost[w]+weightEdge;
                q.push(pair<int,int>(cost[v], v));
            }
        }

    }

}

//Prim & Kruskal

int dj_set(int *color, int s1) {
    if(color[s1]==s1)
        return s1;
    else
        return color[s1] = dj_set(color, color[s1]);
}

void dj_union(int *color,
              int s1, int s2) {
    color[ dj_set(color, s1) ] =
            dj_set(color, s2);
}

int kruskal(Graph *g, Graph *out) {
    initGraph(out, g->V);
    int sum = 0;

    vector< vector<int> > edges(g->E);
    int contEdge = 0;
    for(int u=1; u<=g->V; u++) {
        for(int i=0; i<g->dg[u]; i++) {
            int v = g->adj[u][i];
            int w = g->weight[u][i];
            vector<int> edge(3);
            edge[0] = w;
            edge[1] = u;
            edge[2] = v;
            edges[contEdge++]=edge;
        }
    }
    sort(edges.begin(), edges.end());

    int color[MAXV];
    for(int i=0; i<MAXV; i++)
        color[i] = i;

    for(int i=0; i<g->E; i++) {
#ifdef __DEBUG__
        printf("Aresta de %d para %d com custo %d\n", edges[i][1], edges[i][2], edges[i][0]);
#endif

        int u, v, w;
        w = edges[i][0];
        u = edges[i][1];
        v = edges[i][2];

        if(dj_set(color, u)!=dj_set(color,v)) {
#ifdef __DEBUG__
            printf("Adiciona na solucao!\n");
            printf("ANTES %d %d %d %d %d\n", color[1], color[2], color[3], color[4], color[5]);
#endif
            addEdge(out, u, v, w);
            addEdge(out, v, u, w);
            sum += w;
            /*
            int toChange = color[u];
            for(int k=0; k<MAXV; k++) {
                if(color[k] == toChange)
                    color[k] = color[v];
            }
            */
            dj_union(color, u, v);
#ifdef __DEBUG__
            printf("DPOIS %d %d %d %d %d\n", color[1], color[2], color[3], color[4], color[5]);
#endif
            if(out->E == 2*(g->V-1) ) break;
        }


    }

    return sum;
}
/*
int main() {
    Graph g;
    int v, e;
    scanf("%d%d", &v, &e);

    initGraph(&g, v);

    for(int i=0; i<e; i++) {
        int v, u, w;
        scanf("%d%d%d", &v, &u, &w);
        addEdge(&g, v, u, w);
        addEdge(&g, u, v, w);
    }

    bfs(&g, 1);
    dfs(&g, 1);

    int costs[MAXV];
    dijkstra(&g, 1, costs);

    for(int i=1; i<=g.V; i++) {
        printf("COST[%d] = %d\n", i, costs[i]);
    }

    printf("--- Segunda aula ---\n");

    Graph output;
    printf("Kruskal = %d\n", kruskal(&g, &output) );


}
*/

double distance(vector<int> p1, vector<int> p2){
    double difx = p1[0] - p2[0];
    double dify = p1[1] - p2[1];
    return(sqrt(pow(difx, 2)+pow(dify, 2)));
}

int main() {
    int i, j, n, m, x, y;

    while(scanf("%d", &n)==1){
        vector<vector<int> > pontos(n, vector<int>(2));
        for(i=0; i<n; i++){
            scanf("%d %d", &x, &y);
            pontos[i][0] = x;
            pontos[i][1] = y;
        }
        for(i=0; i<n; i++){
            vector<int> ponto = pontos[i];
            printf("%d, %d\n", ponto[0], ponto[1]);
        }


        Graph g;
        int v = n;

        initGraph(&g, n);
        for(i=0; i<n; i++){
            for(j=i+1; j<n; j++) {
                double d = distance(pontos[i], pontos[j]);
                addEdge(&g, i, j, d);
                addEdge(&g, j, i, d);
            }

        }
        //Printa grafo
        for(int u=0; u<=g.V; u++) {
            for (int i = 0; i < g.dg[u]; i++) {
                int v = g.adj[u][i];
                int w = g.weight[u][i];
                printf ("%d - %d - %d",u,v,w);
            }
        }


//        Graph output;
//        printf("Kruskal = %d\n", kruskal(&g, &output) );
    }
    return 0;
}