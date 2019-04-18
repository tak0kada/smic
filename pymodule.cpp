#include <string>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <array>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cstdlib>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/detail/common.h>
#include <pybind11/embed.h>
namespace py = pybind11;

// #include <eigen3/Eigen/Core>

std::vector<double> split(
        const std::string& data, const char field_delimiter, const int ncol,
        const std::vector<int>& col_idx)
{
    std::vector<double> ret;
    ret.resize(ncol);

    std::string val{};
    int i{0}; // pointer to return vector
    int j{0}; // pointer to current position
    for (const auto& c: data)
    {
        if (j != col_idx[i])
        {
            if (c != field_delimiter)
                continue;
            else
            {
                ++j;
                continue;
            }
        }
        if (c != field_delimiter)
            val += c;
        else
        {
            ret[i] = std::stod(val);
            val = "";
            ++i;
            ++j;
        }
    }
    ret[i] = std::stod(val);
    return ret;
}

struct Mask
{
    Mask(int nrow_file, int ncol_file, const std::vector<int>& skiprows, const std::vector<int>& skipcols)
    {
        row_idx.reserve(nrow_file - skiprows.size());
        col_idx.reserve(ncol_file - skipcols.size());

        // if not blank
        int i{0};
        for (const auto& skiprow: skiprows)
        {
            for (int k = i; k < skiprow; ++k)
                row_idx.push_back(k);
            i = skiprow + 1;
        }
        for (int k = i; k < nrow_file; ++k)
            row_idx.push_back(k);

        // if not blank
        int j{0};
        for (const auto& skipcol: skipcols)
        {
            for (int k = j; k < skipcol; ++k)
                col_idx.push_back(k);
            j = skipcol + 1;
        }
        for (int k = j; k < ncol_file; ++k)
            col_idx.push_back(k);
    }

    std::pair<int, int> operator() (int i, int j)
    {
        return {row_idx[i], col_idx[j]};
    }

    std::vector<int> row_idx;
    std::vector<int> col_idx;
};

struct Config
{
    Config(std::string filename, char field_delimiter, char line_terminator,
           std::vector<int> skiprows = {}, std::vector<int> skipcols = {})
    : filename{filename}, field_delimiter{field_delimiter},
      line_terminator{line_terminator}, ncol_file{0}, nrow_file{0}
    {
        this->skiprows = skiprows;
        this->skipcols = skipcols;

        // open csv file
        std::ifstream ifs{filename};
        if (ifs.fail())
            throw std::runtime_error("Failed to open " + filename);

        // skip BOM
        std::array<char, 3> BOM;
        try
        {
            for (auto& b: BOM)
                ifs.get(b);
        }
        catch (const std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            std::exit(1);
        }
        if (BOM[0] == (char)0xEF && BOM[1] == (char)0xBB && BOM[2] == (char)0xBF)
            ifs.seekg(3, std::ios::beg);
        else
            ifs.seekg(0, std::ios::beg);
        auto beg = ifs.tellg();

        // load all data to memory
        ifs.seekg(0, std::ios::end);
        int filesize = ifs.tellg() - beg;
        auto buf = std::make_unique<char[]>(filesize);
        buf[filesize] = '\0';
        ifs.seekg(beg);
        auto cur = beg;
        ifs.read(buf.get(), filesize);

        // build config
        int ndelim{0};
        int nbreak{0};
        for (int i = 0; i < filesize; ++i)
        {
            if (buf[i] == ',')
            {
                ++ndelim;
            }
            else if (buf[i] == '\n')
            {
                ++nbreak;
                line.push_back({cur, beg + static_cast<std::fstream::pos_type>(i)});
                cur = beg + static_cast<std::fstream::pos_type>(i + 1);
            }
        }
        ncol_file = ndelim / nbreak + 1;
        nrow_file = nbreak;
    }

    std::string filename;
    char field_delimiter;
    char line_terminator;
    int ncol_file;
    int nrow_file;
    std::vector<std::pair<
        std::fstream::pos_type, std::fstream::pos_type>
        > line;
    std::vector<int> skiprows;
    std::vector<int> skipcols;
};

// using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
// RowMatrixXd load_csv(const std::string fname,
// // Eigen::MatrixXd load_csv(const std::string fname,
//          const char delimiter = ',', char line_terminator = '\n',
//          const std::vector<int> skiprows = {}, const std::vector<int> skipcols = {})
// {
//     const Config config{fname, delimiter, line_terminator, skiprows, skipcols};
//     const int nrow = config.nrow_file - skiprows.size();
//     const int ncol = config.ncol_file - skipcols.size();
//     const Mask mask{config.nrow_file, config.ncol_file, config.skiprows, config.skipcols};
    // RowMatrixXd X(nrow, ncol);
    // // Eigen::MatrixXd X(nrow, ncol);
//
// # pragma omp parallel for
//     for (int i = 0; i < nrow; ++i)
//     {
//         int idx = mask.row_idx[i];
//         const auto [beg, end] = config.line[idx];
//         const int length = end - beg;
//
//         std::ifstream ifs(fname);
//         std::string buf(length, ' ');
//         ifs.seekg(beg, std::ios::beg);
//         ifs.read(&buf[0], length);
//
//         const auto data = split(buf, config.field_delimiter, ncol, mask.col_idx);
//         for (int j = 0; j < ncol; ++j)
//         {
//             // X(i, j) = data[j];
//             X.coeffRef(i, j) = data[j];
//         }
//     }
//     return X;
// }

inline py::array_t<double> load_csv(const std::string& fname,
         const char delimiter = ',', char line_terminator = '\n',
         const std::vector<int>& skiprows = {}, const std::vector<int>& skipcols = {})
{
    const Config config{fname, delimiter, line_terminator, skiprows, skipcols};
    const int nrow = config.nrow_file - skiprows.size();
    const int ncol = config.ncol_file - skipcols.size();
    const Mask mask{config.nrow_file, config.ncol_file, config.skiprows, config.skipcols};

    // init np.ndarray
    py::detail::any_container<ssize_t> shape{{nrow, ncol}};
    py::array_t<double> arr{shape};
    auto ref{arr.template mutable_unchecked<2>()};

    bool error_in_loop = false;
    # pragma omp parallel for
    for (int i = 0; i < nrow; ++i)
    {
        int idx = mask.row_idx[i];
        const auto [beg, end] = config.line[idx];
        const int length = end - beg;

        std::ifstream ifs(fname);
        std::string buf(length, ' ');
        ifs.seekg(beg, std::ios::beg);
        ifs.read(&buf[0], length);

        std::unique_ptr<std::vector<double>> data;
        try
        {
            const auto data = split(buf, config.field_delimiter, ncol, mask.col_idx);
            for (int j = 0; j < ncol; ++j)
            {
                ref(i, j) = data[j];
            }
        }
        catch (const std::invalid_argument& e)
        {
            # pragma omp critical
            std::cout << "Error occurred in OpenMP loop, as shown below." << std::endl;
            std::cout << e.what() << std::endl;
            error_in_loop = true;
        }
    }

    if (error_in_loop)
        return {};
    else
        return arr;
}

PYBIND11_MODULE(smic, m)
{
    m.doc() = "C++ based CSV parser binded by pybind11";
    m.def("load_csv", &load_csv, py::return_value_policy::move,
        py::arg("fname"), py::arg("delimiter") = ',', py::arg("line_terminator") = '\n',
        py::arg("skiprows") = std::vector<int>(), py::arg("skipcols") = std::vector<int>());
}
