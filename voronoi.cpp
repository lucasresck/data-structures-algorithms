#include <iostream>
#include <fstream>
#include <queue>

using namespace std;

// TODO:
//  - degenerate case when two sites have the same y coord;

struct point {
    double x, y;
};

typedef point site;

struct circle;

struct arc {
    site s;
    circle *c;
    arc *prev, *next;
};
struct circle{
    double lowest;
    arc *a;
    bool valid;
};

struct compare_events {
    bool operator()(site a, site b) {
        if (a.y < b.y) {
            return true;
        }
        else {
            return false;
        }
    }
    bool operator()(circle *a, circle *b) {
        if (a->lowest < b->lowest) {
            return true;
        }
        else {
            return false;
        }
    }
};

// Event queues
priority_queue<site, vector<site>, compare_events> sites;
priority_queue<circle*, vector<circle*>, compare_events> circles;

int main() {
    ifstream f("points.txt");
    site p;
    while(f >> p.x >> p.y) {
        sites.push(p);
    }
    return 0;
}
