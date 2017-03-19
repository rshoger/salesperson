// vim: set tabstop=4 shiftwidth=4 expandtab :

#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/program_options.hpp>
#include <boost/graph/simple_point.hpp>
#include <boost/graph/metric_tsp_approx.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graphviz.hpp>

using std::cout;
using std::endl;
namespace po = boost::program_options;

/**
 * Process program options using boost::program_options
 */
int process_program_options(int ac, const char** av, po::variables_map *vm)
{
    po::options_description desc("Usage: salesperson [options] [input-file]");

    po::positional_options_description p;

    p.add("input-file", -1);

    desc.add_options()
        ("help,h", "produce help message")
        ("input-file,i", po::value<std::string>(), "supply problem input file")
        ("debug,d", "print debug information")
        ("print,p", "print node information")
    ;

    po::store(po::command_line_parser(ac, av)
            .options(desc).positional(p).run(), *vm);

    po::notify(*vm);

    if (vm->count("help") || ac < 2)
    {
        std::cout << desc << endl;
        return 1;
    }

    return 0;
}

/**
 * Connect the vertices in a graph G. For some vertex u from G create a table
 * such that each cell in the table holds the distance from u to every other
 * vertex in graph G.
 */
template<typename VertexListGraph, typename PointContainer,
         typename WeightMap, typename VertexIndexMap>
void connect_vertex_map(VertexListGraph& g, const PointContainer& points,
        WeightMap wmap,
        VertexIndexMap vmap,
        int /*size*/)
{
    using std::pow;
    using std::sqrt;

    typedef typename boost::graph_traits<VertexListGraph>::edge_descriptor
        Edge;

    typedef typename boost::graph_traits<VertexListGraph>::vertex_iterator
        VIter;

    Edge e;
    bool inserted;

    std::pair<VIter, VIter> verts(vertices(g));
    for (VIter src(verts.first); src != verts.second; src++)
    {
        for (VIter dest(src); dest != verts.second; dest++)
        {
            if (dest != src)
            {
                double weight(sqrt(
                    pow(static_cast<double>(points[vmap[*src]].x -
                        points[vmap[*dest]].x), 2.0) +
                    pow(static_cast<double>(points[vmap[*dest]].y -
                        points[vmap[*src]].y), 2.0)));

                boost::tie(e, inserted) = boost::add_edge(*src, *dest, g);

                wmap[e] = weight;
            }
        }
    }
}

template<typename PositionVector>
void print_position_vector(PositionVector *pos_vec)
{
    using std::vector;
    using boost::simple_point;

    for (vector<simple_point<double> >::iterator itr = pos_vec->begin();
            itr != pos_vec->end(); ++itr++)
        cout << "(" << itr->x << ", " << itr->y << ")" << endl;
}

/**
 * A routine to generate the graph
 */
template<typename PositionVector>
int create_some_graph(PositionVector *position_vec, int n)
{
    using std::vector;

    // Give ourselves a Graph type
    typedef boost::adjacency_matrix<boost::undirectedS, boost::no_property,
        boost::property <boost::edge_weight_t, double> > Graph;

    // Define Vertex type
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

    // Define type VertexContainer to hold elements of type Vertex
    typedef vector<Vertex> VertexContainer;

    // Define type WeightMap to hold a weight map of the graph
    typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

    // Define type VertexMap to hold a vertex map
    typedef boost::property_map<Graph, boost::vertex_index_t>::type VertexMap;

    // A container of vertices
    VertexContainer c;
    // The graph
    Graph g(position_vec->size());

    WeightMap weight_map(get(boost::edge_weight, g));
    VertexMap v_map = get(boost::vertex_index, g);

    connect_vertex_map(g, *position_vec, weight_map, v_map, n);

    print_position_vector(position_vec);

    //lazy short circuit to avoid approx tour
    //return 0;

    boost::metric_tsp_approx_tour(g, std::back_inserter(c));

    double total(0);
    for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr)
    {
        cout << *itr << " ";
        total += boost::lexical_cast<double>(*itr);
    }

    cout << endl << "total = " << total << endl;

    c.clear();

    total = 0;
    boost::metric_tsp_approx_from_vertex(g, *vertices(g).first,
        get(boost::edge_weight, g), get(boost::vertex_index, g),
        boost::tsp_tour_visitor<std::back_insert_iterator<vector<Vertex> > >
        (std::back_inserter(c)));
    for (vector<Vertex>::iterator itr = c.begin(); itr != c.end(); ++itr)
    {
        cout << *itr << " ";
        total += boost::lexical_cast<double>(*itr);
    }
    cout << endl << total << endl;

    return 0;
}

/**
 * Write a file
 */
template<typename Iterator>
int write_a_file(Iterator begin, Iterator end, std::string filename)
{
    //TODO: Ensure that we can take a vector/list of items to be written.
    // The list can be used to supply formatting to a 'dumb' function.

    std::ofstream file(filename.c_str());

    file.exceptions(std::ofstream::badbit);

    try
    {
        for (; begin != end; ++begin)
            file << *begin;

        file << std::endl;
    }
    catch(const std::ofstream::failure& e)
    {
        cout << e.what() << endl;
    }

    file.close();

    return 0;
}

/**
 * Process a file of vertices
 */
int process_problem_file (std::string filename, bool debug)
{
    //TODO: What to do with the first value?

    using std::cout;
    using std::string;
    using std::vector;
    using boost::simple_point;
    using boost::lexical_cast;

    typedef vector<simple_point<double> > PositionVec;

    std::ifstream file(filename.c_str());

    file.exceptions(std::ifstream::badbit);

    PositionVec position_vec;

    int n(0);
    try
    {
        string str;
        while (std::getline(file, str))
        {
            simple_point<double> vertex;

            string x;
            string y;

            boost::char_separator<char> sep(" ");
            boost::tokenizer<boost::char_separator<char> > tokens(str, sep);

            for (auto t = tokens.begin(); t != tokens.end(); ++t)
            {
                auto & v = *t;
                if (t != tokens.begin() && std::next(t) != tokens.end())
                    x = v;
                else if (std::next(t) == tokens.end())
                    y = v;
            }

            if (debug)
            {
                cout << "line: " << str << endl;
                cout << "x: " << x << endl
                     << "y: " << y << endl;
            }

            vertex.x = lexical_cast<double>(x);
            vertex.y = lexical_cast<double>(y);

            position_vec.push_back(vertex);
            n++;
        }
        file.close();
    }
    catch(const std::ifstream::failure& e)
    {
        cout << e.what() << endl;
    }

    if (debug) cout << "Finished processing" << filename << std::endl;

    // TODO: move this outside of this function to main
    create_some_graph(&position_vec, n);
    return 0;
}

int run_simulation (po::variables_map *vm)
{
    // TODO: move graph creation to the 'run' routine

    // Process the problem file supplied in args
    process_problem_file((*vm)["input-file"].as<std::string>(),
            vm->count("debug"));

    // TODO: generate graphviz representation, requires graph reference
    //boost::write_graphviz(std::cout, g);
    return 0;
}

int main (int ac, const char** av)
{
    po::variables_map vm;

    int r = 0;

    // return if we get a failure parsing program options
    if ((r=process_program_options(ac, av, &vm)))
        return r;

    // If an input-file hasn't been provided, display help and exit with error.
    if (!vm.count("input-file"))
    {
        //TODO: un-const this by using char** arg_help = {}
        const char *arg_help[2] = {"0", "--help"};
        process_program_options(2, arg_help, &vm);
        return -1;
    }

    // Run the graph simulation
    run_simulation(&vm);

    // Small example demonstrating how to pass "front" and "back" iterators --
    // begin() and end() respectively -- to a function that writes to a file
    // of a provided filename.
    //
    // A function taking iterators can take iterators of any STL class and
    // and reuse a parsing pattern.
    //
    std::string file_data("something");
    std::string output_filename("some_output.log");
    write_a_file(file_data.begin(), file_data.end(), output_filename);
    return 0;
}
