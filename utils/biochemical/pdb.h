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

/// \file pdb.h
/// Contains class PDB.
/// \author Valerii Sukhorukov

#ifndef UTILS_BIOCHEMICAL_PDB_H
#define UTILS_BIOCHEMICAL_PDB_H

#include <sstream>
#include <algorithm>
#include <string>
#include "../common/misc.h"
#include "../arrays/all.h"
#include "../common/msgr.h"

namespace Utils {
namespace Biochemical {

using namespace Common;
using namespace Arrays;

/// \brief Encapsulates Protein Data Bank formatting and file io.
class Pdb {

public:

	static float scaling;	///< Scaling factor for atom position coordinates.

	A3<float>	pos {huge<float>};				///< Atom position.
	uint		ind {huge<decltype(ind)>};		///< Atom insex.
	int			irecname {};		///< Record name   1 : "ATOM  "; 2 : "HETATM".
	ulong		pdbsegment {};		///< Segment name.
	std::string	name;				///< Pdb name.
	std::string	resname {""};		///< Residue name.
	std::string	chainID {};			///< Chain identifier.
	ulong		resSeq {};			///< Residue sequence number.
	float		occupancy {};		///< Occupancy.
	float		tempFactor {};		///< Temperature factor.
	std::string	element;			///< Name of the element.
	std::string	charge;				///< Electrostatic charge.
	float		vdWRad {};			///< Van Der Waals radius.

	/// \brief Constructor.
	explicit Pdb(const std::string& chainID, const std::string& resname, const std::string& atomname) noexcept;

	/// \brief Constructor.
	explicit Pdb(const uint ind, const std::string& chainID, const uint resSeq, const std::string& resname, const std::string& atomname, const A3<float>& pos) noexcept;

	/// \brief Constructor.
	explicit Pdb(const uint ind, const std::string& chainID, const uint resSeq, const std::string& resname, const std::string& atomname, const A3<float>& pos, const float occupancy) noexcept;

	/// \brief Constructor.
	explicit Pdb(const uint ind, const std::string& chainID, const uint resSeq, const std::string& resname, const std::string& atomname, const A3<float>& pos, const float occupancy, const float tempFactor) noexcept;

	/// \brief Constructor.
	explicit Pdb(const uint ind, const std::string& chainID, const uint resSeq, const std::string& resname, const std::string& atomname, const A3<float>& pos, const float occupancy, const float tempFactor, const std::string& element) noexcept;

	/// \brief Constructor.
	explicit Pdb(const uint ind, const std::string& chainID, const uint resSeq, const std::string& resname, const std::string& atomname, const A3<float>& pos, const std::string& element) noexcept;

	/// \brief Constructor.
	explicit Pdb(const std::string record, const ulong segm) noexcept;

	/// \brief Formated reading of a line record from a pdb file.
	/// \param record Record in the pdb file.
	void format_as_pdb( const std::string& record ) noexcept;

	/// \brief Produce an atom-specific formatted record for a pdb file.
	/// \param msgr Output message processor.
	std::string format_as_pdb( Msgr& msgr ) const;

	/// \brief Produce an atom-specific formatted record for a pdb file.
	/// \param p Atom position coordinates.
	/// \param msgr Output message processor.
	std::string format_as_pdb( const A3<float>& p, Msgr& msgr ) const;

	/// \brief Produce an atom-specific formatted record for a pdb file.
	/// \param p Atom position coordinates.
	/// \param msgr Output message processor.
	std::string format_as_pdb( float* p, Msgr& msgr ) const;

	/// \brief Produce an atom-specific formatted record for a pdb file.
	/// \param a Atom-specific record.
	/// \param msgr Output message processor.
	static std::string format_as_pdb(
		const Pdb& a,
		Msgr& msgr );

	/// \brief Produce an atom-specific formatted record for a pdb file.
	static std::string format_as_pdb(
		const std::string& recn,
		const uint iatom,
		const std::string& atname,
		const std::string& resname,
		const ulong ires,
		const char ichain,
		const A3<float>& pos,
		const float occ,
		const float tempf,
		const std::string& elt,
		const std::string& chge,
		Msgr& msgr );

	/// \brief Determine if the atom is part of an aminoacid.
	bool is_in_aminoacid() const noexcept;

	/// \brief Read from a pdb file.
	/// \param filename Name of the file to read from.
	/// \param a Atom-specific records.
	/// \param other Non-atom records.
	/// \param msgr Output message processor.
	static void read( const std::string& filename,
					  std::vector<Pdb>& a,
		  		 	  vec2<std::string>& other,
		  		 	  Msgr& msgr );

	/// \brief Write to a pdb file.
	/// \param filename Name of the file to write to.
	/// \param a Atom-specific records.
	/// \param other Non-atom records.
	/// \param msgr Output message processor.
	static void write( const std::string& filename,
				  	   const std::vector<Pdb>& a,
				  	   const vec2<std::string>& other,
				  	   Msgr& msgr );
	   
private:


	/// \brief Convert record name string to index.
	/// \param s "HETATM" or "ATOM  ".
	int irecn(const std::string& s) const noexcept;

	/// \brief Get record name.
	/// \return "HETATM" or "ATOM  ".
	std::string recname() const noexcept; 

	/// \brief Convert record name index to string.
	/// \return "HETATM" or "ATOM  ".
	static std::string irecname2str(const int irecname) noexcept;
	
	/// \brief Determine van Der Waals radius based on the atom data.
	float set_vdW() const noexcept;

};

}	// namespace Biochemical
}	// namespace Utils

#endif // UTILS_BIOCHEMICAL_PDB_H








