/* =============================================================================

 Copyright (C) 2009-2025 Valerii Sukhorukov. All Rights Reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

================================================================================
*/

/**
 * \file dx_io.h
 * \brief Functions reading and writing volumetric data using .DX file format.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_VOLUMETRIC_DX_IO
#define UTILS_VOLUMETRIC_DX_IO

#include <concepts>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "../arrays/all.h"
#include "../common/misc.h"
#include "../constants.h"
#include "../msgr.h"

/// Volume data io.
namespace utils::volumetric {

using namespace constants;

/**
 * \brief Writes 3d matrix to a file using .DX format.
 * \tparam K Type of the original data values.
 * \tparam L Floating type used for recording.
 * \tparam T Floating type used in dimensional description.
 * \param filename Name of the file.
 * \param v 3d matrix containing data to save.
 * \param origin Position of the coordinate system origin.
 * \param delta Step sizes in each dimension.
 * \param units Length units used.
 * \param msgr Output message processor.
 */
template<std::floating_point K,
         std::floating_point L,
         arithmetic T>
void save_as_DX(
    const std::string& filename,
    const vec3<K>& v,
    arrays::A3<T> origin,        // by value
    arrays::A3<T> delta,         // by value
    const std::string& units,
    Msgr& msgr
)
{
    if (units == "nm") {
        delta *= ten<T>;
        origin *= ten<T>;
    }
    else if (units != "A")
        msgr.exit("Error while saving volume map: length units not supported");
    std::ofstream fout {filename};
    if (fout.fail())
        msgr.exit("Error while saving volume map: Unable to open file: " +
                  filename);

    msgr.print("Started writing to " + filename);

    std::size_t npoints = v.size() * v[0].size() * v[0][0].size();

    fout << "object 1 class gridpositions counts "
         << v.size() << " " << v[0].size() << " " << v[0][0].size()
         << std::endl;

    fout << std::string("origin ")
         << origin[0] << " " << origin[1] << " " << origin[2]
         << std::endl;

    fout << "delta " << delta[0] << " " << 0        << " " << 0        << std::endl;
    fout << "delta " << 0        << " " << delta[1] << " " << 0        << std::endl;
    fout << "delta " << 0        << " " << 0        << " " << delta[2] << std::endl;

    fout << "object 2 class gridconnections counts "
         << v.size() << " " <<v[0].size() << " " << v[0][0].size()
         << std::endl;
    fout << "object 3 class array type double rank 0 items "
         << npoints << " data follows"
         << std::endl;

    std::vector<L> templine(npoints);
    std::size_t count {};
    for (std::size_t i=0; i<v.size(); i++)
        for (std::size_t j=0; j<v[i].size(); j++)
            for (std::size_t k=0; k<v[i][j].size(); k++)
                templine[count++] = v[i][j][k];
    std::size_t nlines = std::size_t(npoints/3);
    std::size_t lastlinenum = npoints % 3;
    for (size_t i=0; i<nlines; i++)
        fout << STR(templine[3*i]) + " " +
                STR(templine[3*i+1]) + " " +
                STR(templine[3*i+2]) << std::endl;
    if (     lastlinenum == 1) fout << STR(templine[npoints-1]) << std::endl;
    else if (lastlinenum == 2) fout << STR(templine[npoints-2]) + " " +
                                       STR(templine[npoints-1]) << std::endl;
    fout << std::endl
         << "object \"distance (protein) [A]\" class field"
         << std::endl
         << std::endl;

    msgr.print("Finished writing to " + filename);
}

/**
* \brief Reads 3d volumetric data from a .DX format file into a matrix.
* \tparam K Type of the output data values.
* \tparam T Floating type used in the .DX records.
* \param filename Name of the file.
* \param[out] v 3d matrix containing data to save.
* \param[out] origin Position of the coordinate system origin.
* \param[out] delta Step sizes in each dimension.
* \param[in] bin Ddata dimensions.
* \param[in] units Length units used.
* \param msgr Output message processor.
*/
template<arithmetic K,
         std::floating_point T>
void read_as_DX(
    const std::string& filename,
    vec3<K>& v,
    arrays::A3<T>& origin,
    arrays::A3<T>& delta,
    arrays::A3<szt>& ms,
    const std::vector<uint>& bin,
    const std::string& units,
    Msgr &msgr
)
{
    if (bin[0] != 1 ||
        bin[1] != 1 ||
        bin[2] != 1)
        msgr.exit("Error: readDX with bin != 1 is not implemented for " +
                  filename);

    std::ifstream fin {filename, std::ios::in};
    if (fin.fail())
        msgr.exit("Error while reading volume map: Unable to open file at " +
                  filename);

    msgr.print("Started reading " + filename);

    std::string line {};
    szt xdim = 0;
    szt ydim = 0;
    szt zdim = 0;
    // % object 1 class gridpositions counts 30 30 30:
    std::getline(fin, line);
    sscanf(line.c_str(),
           "%*s %*d %*s %*s %*s %ld %ld %ld",
           &xdim, &ydim, &zdim);
    ms[0] = xdim;
    ms[1] = ydim;
    ms[2] = zdim;
    v = common::Vec3::make<K>(xdim, ydim, zdim);
    auto npoints {xdim * ydim * zdim};
    szt numlines;
    auto lastlinenumelements = npoints % 3;
    if (lastlinenumelements == 0) {
        lastlinenumelements = 3;
        numlines = npoints / 3;
    }
    else
        numlines = static_cast<szt>(npoints / 3) + 1;

    // origin:
    std::getline(fin, line);
    sscanf(line.c_str(), "%*s %f %f %f", &origin[0], &origin[1], &origin[2]);

    // delta 1 0 0:
    std::getline(fin, line);
    sscanf(line.c_str(), "%*s %f %*f %*f", &delta[0]);

    // delta 0 1 0:
    std::getline(fin, line);
    sscanf(line.c_str(), "%*s %*f %f %*f", &delta[1]);

    // delta 0 0 1:
    std::getline(fin, line);
    sscanf(line.c_str(), "%*s %*f %*f %f", &delta[2]);

    // Ignore: object 2 class grid connections counts 247 317 242:
    std::getline(fin, line);

    // Ignore: object 3 class array type double rank 0 items 18948358 data follows:
    std::getline(fin, line);

    // Convert to float for writing:
    std::vector<float> templine(npoints);

    szt count {};
    for (szt n=0; n<numlines-1; n++) {
        std::getline(fin, line);
        sscanf(line.c_str(),
               "%f %f %f",
               &templine[count], &templine[count+1], &templine[count+2]);
        count += 3;
    }

    // The last line may be shorter than 3 elements, so it is treated separately:
    std::getline(fin, line);
    if (     lastlinenumelements == 1) {
        sscanf(line.c_str(),
               "%f",
               &templine[count]);
        count++;
    }
    else if (lastlinenumelements == 2) {
        sscanf(line.c_str(),
               "%f %f",
               &templine[count], &templine[count+1]);
        count += 2;
    }
    else {
        sscanf(line.c_str(),
               "%f %f %f",
               &templine[count], &templine[count+1], &templine[count+2]);
        count += 3;
    }

    count = 0;
    for (szt i=0; i<v.size(); i++)
        for (szt j=0; j<v[i].size(); j++)
            for (szt k=0; k<v[i][j].size(); k++)
                v[i][j][k] = static_cast<K>(templine[count++]);
    if (units == "nm") {
        delta /= ten<T>;
        origin /= ten<T>;
    }
    else if (units != "A")
        msgr.exit("Error while reading volume map: length units not supported ");

    msgr.print("Finished reading " + filename);
}

}  // namespace utils::volumetric

#endif  // UTILS_VOLUMETRIC_DX_IO
