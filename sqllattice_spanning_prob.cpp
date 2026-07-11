#include <bits/stdc++.h>
using namespace std;

const double EPS_COORD = 1e-9;

struct ROW {
    int id;
    vector<tuple<double, double, double, double>> points;
    vector<int> up_neighbors;    
    vector<int> down_neighbors;
    vector<int> left_neighbors;
    vector<int> right_neighbors;
    double area;
};

static thread_local std::mt19937_64 rng((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());

int findroot(int x, vector<int>& parent) {
    if (parent[x] < 0) return x;
    return parent[x] = findroot(parent[x], parent);
}

void unite(int a, int b, vector<int>& parent, vector<int>& boundary_mask) {
    a = findroot(a, parent);
    b = findroot(b, parent);
    if (a == b) return;

    if (parent[a] < parent[b]) {  
        parent[a] += parent[b];
        parent[b] = a;
        boundary_mask[a] |= boundary_mask[b];
    } else {
        parent[b] += parent[a];
        parent[a] = b;
        boundary_mask[b] |= boundary_mask[a];
    }
}

int percolator(int L_side, const vector<ROW>& rows, bool save_files) {
    int total_cells = L_side * L_side;

    ofstream tout;
    ofstream chosenFile;
    if (save_files) {
        tout.open("percolation_n.txt");
        chosenFile.open("chosen_cells.txt");
    }

    vector<int> parent(total_cells, 0);
    vector<int> boundary_mask(total_cells, 0);
    vector<bool> occupied(total_cells, false);
    vector<vector<int>> nn(total_cells);

    const int TOUCH_TOP    = 1 << 0;
    const int TOUCH_BOTTOM = 1 << 1;
    const int TOUCH_LEFT   = 1 << 2;
    const int TOUCH_RIGHT  = 1 << 3;

    for (int i = 0; i < total_cells; ++i) {
        const ROW& cell = rows[i];
        nn[i].reserve(cell.up_neighbors.size() + cell.down_neighbors.size() + cell.left_neighbors.size() + cell.right_neighbors.size());
        for (int n : cell.up_neighbors)    nn[i].push_back(n);
        for (int n : cell.down_neighbors)  nn[i].push_back(n);
        for (int n : cell.left_neighbors)  nn[i].push_back(n);
        for (int n : cell.right_neighbors) nn[i].push_back(n);

        if (save_files) {
            tout << i << "\t";
            for (size_t j = 0; j < nn[i].size(); ++j) {
                tout << nn[i][j];
                if (j != nn[i].size() - 1) tout << ",";
            }
            tout << "\n";
        }
    }

    vector<int> order(total_cells);
    iota(order.begin(), order.end(), 0);
    shuffle(order.begin(), order.end(), rng);

    int spanning_cluster_size = 0;

    for (int t = 0; t < total_cells; ++t) {
        int s = order[t];
        occupied[s] = true;
        parent[s] = -1;

        if (save_files) {
            chosenFile << rows[s].id << endl;
            for (const auto& point : rows[s].points) {
                double y_min, y_max, x_min, x_max;
                tie(y_min, y_max, x_min, x_max) = point;
                chosenFile << x_min << " " << y_min << "\n" << x_max << " " << y_min << "\n"
                           << x_max << " " << y_max << "\n" << x_min << " " << y_max << "\n\n";
            }
        }

        auto& p = rows[s].points[0];
        double y_min = get<0>(p);
        double y_max = get<1>(p);
        double x_min = get<2>(p);
        double x_max = get<3>(p);

        if (fabs(y_min - 0.0) < EPS_COORD)          boundary_mask[s] |= TOUCH_BOTTOM;
        if (fabs(y_max - (double)L_side) < EPS_COORD) boundary_mask[s] |= TOUCH_TOP;
        if (fabs(x_min - 0.0) < EPS_COORD)          boundary_mask[s] |= TOUCH_LEFT;
        if (fabs(x_max - (double)L_side) < EPS_COORD) boundary_mask[s] |= TOUCH_RIGHT;

        for (int nb : nn[s]) {
            if (occupied[nb]) {
                unite(s, nb, parent, boundary_mask);
            }
        }

        int root = findroot(s, parent);
        int mask = boundary_mask[root];

        if (((mask & TOUCH_LEFT) && (mask & TOUCH_RIGHT)) || 
            ((mask & TOUCH_TOP) && (mask & TOUCH_BOTTOM))) {
            spanning_cluster_size = t + 1;
            break;
        }
    }

    if (save_files) {
        tout.close();
        chosenFile.close();
    }

    return spanning_cluster_size;
}

vector<ROW> generate_square_lattice(int L_side) {
    int total_cells = L_side * L_side;
    vector<ROW> rows(total_cells);

    for (int r = 0; r < L_side; ++r) {
        for (int c = 0; c < L_side; ++c) {
            int id = r * L_side + c;
            
            double y_min = (double)r;
            double y_max = (double)(r + 1);
            double x_min = (double)c;
            double x_max = (double)(c + 1);

            rows[id].id = id;
            rows[id].points = {make_tuple(y_min, y_max, x_min, x_max)};
            rows[id].area = 1.0;

            if (r > 0)           rows[id].down_neighbors.push_back((r - 1) * L_side + c);
            if (r < L_side - 1)  rows[id].up_neighbors.push_back((r + 1) * L_side + c);
            if (c > 0)           rows[id].left_neighbors.push_back(r * L_side + (c - 1));
            if (c < L_side - 1)  rows[id].right_neighbors.push_back(r * L_side + (c + 1));
        }
    }
    return rows;
}

int main() {
    int L_side = 71; 
    int total_cells = L_side * L_side;  
    int iterations = 10000;   

    vector<double> percolation_stats(total_cells + 2, 0.0);

    for (int trial = 1; trial <= iterations; ++trial) {
        vector<ROW> rows = generate_square_lattice(L_side);  
        
        bool save_files = (trial == iterations); 
        int cells_to_span = percolator(L_side, rows, save_files);  

        if (save_files) {
            ofstream latticeFile("lattice.txt");
            for (const auto& row : rows) {
                latticeFile << row.id << "\n";
                for (const auto& point : row.points) {
                    double y_min, y_max, x_min, x_max;
                    tie(y_min, y_max, x_min, x_max) = point;
                    latticeFile << x_min << " " << y_min << "\n" << x_max << " " << y_min << "\n"
                                << x_max << " " << y_max << "\n" << x_min << " " << y_max << "\n"
                                << x_min << " " << y_min << "\n";
                }
                latticeFile << "\n";
            }
            latticeFile.close();
        }

        if (cells_to_span > 0 && cells_to_span <= total_cells) {
            for (int g = cells_to_span; g <= total_cells; ++g) {
                percolation_stats[g] += 1.0;
            }
        }

        if (trial % max(1, iterations / 10) == 0) {
            cout << "Completed " << trial << " / " << iterations << endl;
        }
    }

    ofstream fout("percolation_results8000.txt");
    fout << "FractionOccupied Frequency" << endl;
    for (int i = 1; i <= total_cells; ++i) {
        fout << (double(i) / total_cells) << " " << (percolation_stats[i] / iterations) << "\n";
    }
    fout.close();

    return 0;
}
