/* ==============================================================================

 Copyright (C) 2009-2019, Valerii Sukhorukov, <vsukhorukov@yahoo.com>

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

============================================================================== */

/// \file volume_io.h
/// \brief Functions reading and writing volumetric data using files using .DX format.
/// \author Valerii Sukhorukov

#ifndef UTILS_VOLUME_IO
#define UTILS_VOLUME_IO

#include "misc.h"
#include "msgr.h"

/// Library-wide.
namespace Utils {
/// General stuff.
namespace Common {

/// \brief Writes 3d matrix to a file using .DX format.
/// \tparam K Type of the original data values.
/// \tparam L Floating type used for recording.
/// \tparam T Floating type used in dimensional description.
/// \param filename Name of the file.
/// \param v 3d matrix containing data to save.
/// \param origin Position of the coordinate system origin.
/// \param delta Step sizes in each dimension.
/// \param units Length units used.
/// \param msgr Output message processor.
template <typename K, typename L, typename T>
static void save_as_DX(
	const std::string& filename,
	const vec3<K>& v,
	A3<T> origin,		// by value
	A3<T> delta,		// by value
	const std::string& units,
	Msgr& msgr )
{
	if (units == "nm") {
		delta *= ten<T>;
		origin *= ten<T>;
	}
	else if (units != "A")
		msgr.exit("Error while saving volume map: length units not supported");
	std::ofstream fout {filename};
	if (!fout.is_open())
		msgr.exit("Error while saving volume map: Unable to open file: "+filename);
	std::cout << "Started writing to "+filename << std::endl;
	std::size_t npoints = v.size() * v[0].size() * v[0][0].size();

	fout << std::string( "object 1 class gridpositions counts ") + STR(v.size()) + " "
																 + STR(v[0].size()) + " "
																 + STR(v[0][0].size()) << std::endl;
	fout << std::string("origin ") + STR(origin[0]) + " " + STR(origin[1]) + " " + STR(origin[2]) << std::endl;

	fout << std::string("delta ") + STR(delta[0]) + " " + STR(0)        + " " + STR(0) << std::endl;
	fout << std::string("delta ") + STR(0)        + " " + STR(delta[1]) + " " + STR(0) << std::endl;
	fout << std::string("delta ") + STR(0)        + " " + STR(0)        + " " + STR(delta[2]) << std::endl;

	fout << std::string("object 2 class gridconnections counts ") + STR(v.size()) + " "
																  + STR(v[0].size()) + " "
																  + STR(v[0][0].size()) << std::endl;
	fout << std::string("object 3 class array type double rank 0 items ") + STR(npoints) + " data follows" << std::endl;

	std::vector<L> templine(npoints);
	std::size_t count {};
	for (std::size_t i=0; i<v.size(); i++)
		for (std::size_t j=0; j<v[i].size(); j++)
			for (std::size_t k=0; k<v[i][j].size(); k++)
				templine[count++] = v[i][j][k];				// ( v[i][j][k] ) ? 0 : 1;
	std::size_t nlines = std::size_t(npoints/3);
	std::size_t lastlinenum = npoints % 3;
	for (size_t i=0; i<nlines; i++)
		fout << STR(templine[3*i])+" "+
				STR(templine[3*i+1])+" "+
				STR(templine[3*i+2]) << std::endl;
	if (     lastlinenum == 1) fout << STR(templine[npoints-1]) << std::endl;
	else if (lastlinenum == 2) fout << STR(templine[npoints-2]) + " " +
									   STR(templine[npoints-1]) << std::endl;
	fout << std::endl;
	fout << std::string("object \"distance (protein) [A]\" class field") << std::endl;
	fout << std::endl;
	std::cout << "Finished writing to "+filename << std::endl;
}

/// \brief Reads 3d volumetric data from a .DX format file into a matrix.
/// \tparam K Type of the output data values.
/// \tparam T Floating type used in the .DX records.
/// \param filename Name of the file.
/// \param[out] v 3d matrix containing data to save.
/// \param[out] origin Position of the coordinate system origin.
/// \param[out] delta Step sizes in each dimension.
/// \param[in] bin Ddata dimensions.
/// \param[in] units Length units used.
/// \param msgr Output message processor.
template <typename K, typename T>
static void readDX( const std::string& filename,
					vec3<K>& v,
					A3<T>& origin,
					A3<T>& delta,
					A3<szt>& ms,
					const std::vector<uint>& bin,
					const std::string& units,
					Msgr &msgr )
{
	if (bin[0] != 1 ||
		bin[1] != 1 ||
		bin[2] != 1 )
		msgr.exit("Error: readDX with bin != 1 is not implemented for "+filename);
	std::ifstream fin {filename, std::ios::in};
	if (!fin.is_open())
		msgr.exit("Error while reading volume map: Unable to open file at "+filename);

	msgr.print{"Started reading "+filename);
	std::string line;
	szt xdim, ydim, zdim;
	std::getline(fin, line);
	sscanf(line.c_str(), "%*s %*d %*s %*s %*s %ld %ld %ld", &xdim, &ydim, &zdim);	// % object 1 class gridpositions counts 30 30 30
	ms[0] = xdim;
	ms[1] = ydim;
	ms[2] = zdim;
	v = make_vec3<K>(xdim, ydim, zdim);
	auto npoints {xdim * ydim * zdim};
	szt numlines;
	auto lastlinenumelements = npoints % 3;
	if (lastlinenumelements == 0) {
		lastlinenumelements = 3;
		numlines = npoints / 3;
	}
	else
		numlines = size_t(npoints/3) + 1;
	std::getline(fin, line);
	sscanf(line.c_str(), "%*s %f %f %f", &origin[0], &origin[1], &origin[2]);	// origin
	std::getline(fin, line);
	sscanf(line.c_str(), "%*s %f %*f %*f", &delta[0]);							// delta 1 0 0
	std::getline(fin, line);
	sscanf(line.c_str(), "%*s %*f %f %*f", &delta[1]);							// delta 0 1 0
	std::getline(fin, line);
	sscanf(line.c_str(), "%*s %*f %*f %f", &delta[2]);							// delta 0 0 1
	std::getline(fin, line);			// ignore:   object 2 class gridconnections counts 247 317 242
	std::getline(fin, line);			// ignore:   object 3 class array type double rank 0 items 18948358 data follows
	std::vector<float> templine(npoints);	// convert to float for writing

	szt count {};
	for (szt n=0; n<numlines-1; n++) {
		std::getline(fin, line);
		sscanf(line.c_str(), "%f %f %f", &templine[count], &templine[count+1], &templine[count+2]);
		count += 3;
	}
	std::getline(fin, line);		// the last line may be shorter than 3 elements, so it is treated separately
	if (     lastlinenumelements == 1)
		sscanf(line.c_str(), "%f", &templine[count++]);
	else if (lastlinenumelements == 2 ) {
		sscanf(line.c_str(), "%f %f", &templine[count], &templine[count+1]);
		count += 2;
	}
	else {
		sscanf(line.c_str(), "%f %f %f", &templine[count], &templine[count+1], &templine[count+2]);
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
		msgr.exit("Error while saving volume map: length units not supported  ");

	msgr.print("Finished reading " + filename);
}

} 	// namespace Common
}	// namespace Utils

#endif // UTILS_VOLUME_IO
