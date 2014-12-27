#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

double curr_x = 0;
double eps = 10e-7;

struct Segment;

struct Point {
    double x;
    double y;

    string _orient;
    Segment *_s1 = nullptr;
    Segment *_s2 = nullptr;

    Point(double x, double y): x(x), y(y) {}
    Point(): Point(0, 0) {}

    friend ostream &operator<<(ostream &os, const Point &p);
};

bool operator==(const Point &self, const Point &other) {
    return self.x == other.x && self.y == other.y;
}

bool operator<(const Point &self, const Point &other) {
    if (self.x < other.x) return true;
    if (self.x > other.x) return false;
    if (self.y < other.y) return true;
    if (self.y > other.y) return false;
    return false;
}

bool operator>(const Point &self, const Point &other) {
    if (self.x > other.x) return true;
    if (self.x < other.x) return false;
    if (self.y > other.y) return true;
    if (self.y < other.y) return false;
    return false;
}

Point operator-(const Point &self, const Point &other) {
    return Point(self.x - other.x, self.y - other.y);
}

ostream &operator<<(ostream &os, const Point &p) {
    os << "(" << p.x << "," << p.y << ")";
    return os;
}

struct Segment {
    Point p1;
    Point p2;
    int ind;

    Segment(Point &&p1, Point &&p2, int ind = 0): p1(p1), p2(p2),
        ind(ind) {}

    Segment(): Segment(Point(), Point(), 0) {}

    double compute_value() {
        Point diff(p2 - p1);
        return p1.y + (curr_x - p1.x) * diff.y / diff.x;
    }

    friend ostream &operator<<(ostream &os, const Segment &s);
};

ostream &operator<<(ostream &os, const Segment &s) {
    os << s.p1 << " -> " << s.p2 << " " << s.ind << endl;
    return os;
}

bool operator==(const Segment &self, const Segment &other) {
    Point diff(self.p2 - self.p1);
    double val1 = self.p1.y + (curr_x - self.p1.x) * diff.y / diff.x;
    diff = other.p2 - other.p1;
    double val2 = other.p1.y + (curr_x - other.p1.x) * diff.y / diff.x;
    if (abs(val1 - val2) < eps)  return true;
    return false;  
}

bool operator<(const Segment &self, const Segment &other) {
    Point diff(self.p2 - self.p1);
    double val1 = self.p1.y + (curr_x - self.p1.x) * diff.y / diff.x;
    diff = (other.p2 - other.p1);
    double val2 = other.p1.y + (curr_x - other.p1.x) * diff.y / diff.x;
    if (val1 < val2) return true;
    if (val1 > val2) return false;
    return false;
}

map<pair<Point, Point>, double> keydict;
set<pair<double, double>> already_enqued;
map<pair<double, double>, pair<Segment, Segment>> intersection_segment_map;

bool intersect(const Segment *s1, const Segment *s2) {
    if (s1 == nullptr || s2 == nullptr) {
        return false;
    }
    const Point &A = s1->p1;
    const Point &B = s1->p2;
    const Point &C = s2->p1;
    const Point &D = s2->p2;
    auto ccw = [](const Point &A, const Point &B, const Point &C) {
        return (C.y - A.y)*(B.x - A.x) > (B.y - A.y)*(C.x - A.x);
    };
    return ccw(A, C, D) != ccw(B, C, D) && ccw(A, B, C) != ccw(A, B, D);
} 

void debug_sweep_line(map<Segment, Segment> &SL) {
    cout << "Sweep Line:" << endl;
    vector<Segment> l;
    for (auto &node: SL) {
        l.push_back(node.first);
    }
    for (auto it = l.rbegin(); it != l.rend(); ++it) {
        Segment s = *it;
        cout << s.p1.x << " " << s.p1.y << " " << s.p2.x << " " << s.p2.y 
             << " " << s.ind << endl;
    }
    cout << "---------------------" << endl;
}

const Segment *succ(map<Segment, Segment> &SL, Segment &s) {
    auto e = SL.end();
    --e;
    auto it = SL.find(s);
    if (it->first == e->first) {
        return nullptr;
    }
    ++it;
    return &it->first;
}

const Segment *pred(map<Segment, Segment> &SL, Segment &s) {
    auto b = SL.begin();
    auto it = SL.find(s);
    if (it->first == b->first) {
        return nullptr;
    }
    --it;
    return &it->first;
}

void handle_intersect(const Segment *s1, const Segment *s2, 
                      map<Segment, Segment> &SL, map<Point, Segment> &p_s, 
                      priority_queue<Point, vector<Point>, greater<Point>> &EQ);

void handle_intersect_event(const Segment *s1, const Segment *s2,
                            map<Segment, Segment> &SL, map<Point, Segment> &p_s,
                            priority_queue<Point, vector<Point>, greater<Point>> &EQ) {
    double temp = curr_x;
    if (s1 == nullptr || s2 == nullptr) {
        return;
    }
    Segment seg1 = const_cast<Segment&>(*s1);
    Segment seg2 = const_cast<Segment&>(*s2);
    auto seg1_pair = make_pair(seg1.p1, seg1.p2);
    auto seg2_pair = make_pair(seg2.p1, seg2.p2);
    curr_x = keydict[seg1_pair];
    auto it = SL.find(seg1);
    if (it == SL.end()) {
        curr_x = temp;
        auto it2 = SL.find(seg1);
        if (it2 == SL.end()) {
        } else {
            SL.erase(it2);
        }
    } else {
        SL.erase(it);
    }
    curr_x = keydict[seg2_pair];
    SL.erase(seg2);
    curr_x = temp + 0.001;
    SL.emplace(seg1, seg1);
    keydict[seg1_pair] = curr_x;
    SL.emplace(seg2, seg2);
    keydict[seg2_pair] = curr_x;
    Segment mini, maxi;
    if (seg1 < seg2) {
        mini = seg1;
        maxi = seg2;
    } else {
        mini = seg2;
        maxi = seg1;
    }
    const Segment *below_mini = pred(SL, mini);
    const Segment *above_maxi = succ(SL, maxi);
    if (intersect(below_mini, &mini)) {
        handle_intersect(below_mini, &mini, SL, p_s, EQ);
    }
    if (intersect(above_maxi, &maxi)) {
        handle_intersect(above_maxi, &maxi, SL, p_s, EQ);
    }
    curr_x = temp;
}

void handle_intersect(const Segment *s1, const Segment *s2,
                      map<Segment, Segment> &SL, map<Point, Segment> &p_s,
                      priority_queue<Point, vector<Point>, greater<Point>> &EQ) {
    double x1 = s1->p1.x;
    double y1 = s1->p1.y;
    double x2 = s1->p2.x;
    double y2 = s1->p2.y;
    double x3 = s2->p1.x;
    double y3 = s2->p1.y;
    double x4 = s2->p2.x;
    double y4 = s2->p2.y;
    double X = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) /
               ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
    double Y = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) /
               ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
    Point point(X, Y);
    auto p = make_pair(X, Y);
    if (already_enqued.find(p) == already_enqued.end()) {
        already_enqued.insert(p);
        auto sp = make_pair(*s1, *s2);
        intersection_segment_map.emplace(p, sp);
        point._s1 = new Segment(*const_cast<Segment*>(s1));
        point._s2 = new Segment(*const_cast<Segment*>(s2));
        point._orient = "intersect";
        EQ.push(point);
    }
}

bool shamos_hoey(vector<Segment> &segments, double lag=0.01) {
    map<Segment, Segment> SL;
    priority_queue<Point, vector<Point>, greater<Point>> EQ;
    cout << "Total number of segments: " << segments.size() << endl;
    map<Point, Segment> p_s;
    for (auto &s: segments) {
        p_s.emplace(s.p1, s);
        p_s.emplace(s.p2, s);
    }
    for (auto &s: segments) {
        s.p1._orient = "left";
        s.p2._orient = "right";
        EQ.emplace(s.p1);
        EQ.emplace(s.p2);
    }
    bool destroy = false;
    while (EQ.size() > 0) {
        Point E = EQ.top();
        string orient = E._orient;
        if (orient == "left") {
            Segment segE = p_s[E];
            curr_x = segE.p1.x;
            SL.emplace(segE, segE);
            keydict[make_pair(segE.p1, segE.p2)] = curr_x;
            const Segment *segA = succ(SL, segE);
            const Segment *segB = pred(SL, segE);
            if (intersect(&segE, segA)) {
                handle_intersect(&segE, segA, SL, p_s, EQ);
                if (destroy) return true;
            }
            if (intersect(&segE, segB)) {
                handle_intersect(&segE, segB, SL, p_s, EQ);
                if (destroy) return true;
            }
        } else if (orient == "right") {
            Segment segE = p_s[E];
            const Segment *segA = succ(SL, segE);
            const Segment *segB = pred(SL, segE);
            curr_x = keydict[make_pair(segE.p1, segE.p2)];
            SL.erase(segE);
            curr_x = segE.p2.x;
            if (intersect(segA, segB)) {
                handle_intersect(segA, segB, SL, p_s, EQ);
                if (destroy) return true;
            }
        } else if (orient == "intersect") {
            curr_x = E.x;
            Segment *s1 = E._s1;
            Segment *s2 = E._s2;
            handle_intersect_event(s1, s2, SL, p_s, EQ);
        }
        EQ.pop();
    }
    return false;
}

void summary() {
    int total = 0;
    for (auto &intersection: already_enqued) {
        ++total;
        auto s = intersection_segment_map[intersection];
        auto s1 = s.first;
        auto s2 = s.second;
        cout << "----------------------------" << endl;
        cout << "Intersection: " << " (" << intersection.first << "," << 
                intersection.second << ")" << endl;
        cout << "Segment1: " << " (" << s1.p1.x << "," << s1.p1.y << ")" << 
            " => " << "(" << s1.p2.x << "," << s1.p2.y << ")" << endl;
        cout << "Segment2: " << " (" << s2.p1.x << "," << s2.p1.y << ")" << 
            " => " << "(" << s2.p2.x << "," << s2.p2.y << ")" << endl;
        cout << "Total: " << total << endl;
    }
}

int main(int argc, char *argv[])
{
    vector<Segment> segments {
        Segment(Point(252.8422544759612, 267.0866807045632),
                Point(306.56700684941103, 149.86014154341402), 0),
        Segment(Point(91.31754664594935, 91.31754664594935),
                Point(327.44042622823946, 64.08873524224386), 1),
        Segment(Point(6.243151567232763, 6.243151567232763),
                Point(84.3195455983842, 34.0661824429703), 2),
        Segment(Point(311.840265890957, 165.3625530190342),
                Point(341.5941181926875, 140.29848154789053), 3),
        Segment(Point(73.10505297848016, 270.4695928644747),
                Point(201.87440935815667, 201.87440935815667), 4),
        Segment(Point(48.24348894060297, 87.24516000564093),
                Point(48.80876112959429, 387.5018218802504), 5),
        Segment(Point(182.41039593940553, 196.4224103209281),
                Point(327.02971692603853, 272.4898099868745), 6),
        Segment(Point(54.847240442794075, 311.0406370779017),
                Point(247.29918986634468, 161.31338165137024), 7),
        Segment(Point(104.66634908254039, 94.07145379266844),
                Point(109.60847065923866, 109.60847065923866), 8),
        Segment(Point(48.454472183760174, 48.454472183760174),
                Point(117.87458829812141, 288.99602998668644), 9),
        Segment(Point(108.31439308252064, 222.5956000814636),
                Point(195.59909368691248, 190.5910930136407), 10),
        Segment(Point(63.14980626808637, 198.71552141870205),
                Point(380.37347160181554, 377.06527181419386), 11),
        Segment(Point(108.58324630049717, 289.4927081947947),
                Point(398.9906017217017, 192.37254161116445), 12),
        Segment(Point(257.7006344089656, 283.3226278498322),
                Point(334.6937711655362, 113.94721821399449), 13),
        Segment(Point(35.85236105711878, 241.8903816851933),
                Point(216.71648917784938, 113.98872400541632), 14)
    };
    cout << shamos_hoey(segments) << endl;
    summary();
    return 0;
}
