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

============================================================================== */

/**
 * \file pdb.h
 * \brief Contains class PDB.
 * \author Valerii Sukhorukov
 */


#ifndef UTILS_BIOCHEMICAL_PDB_H
#define UTILS_BIOCHEMICAL_PDB_H

#include <algorithm>
#include <filesystem>
#include <sstream>
#include <string>

#include "../arrays/all.h"
#include "../common/misc.h"
#include "../constants.h"
#include "../msgr.h"

/// Biochemistry-related stuff.
namespace utils::biochemical {

/// Encapsulates Protein Data Bank formatting and file io.
class Pdb {

public:

    using A3f = utils::arrays::A3<float>;

    static float scaling;    ///< Scaling factor for atom position coordinates.

    A3f           pos {undefined<float>};         ///< Atom position.
    unsigned int  ind {undefined<unsigned int>};  ///< Atom insex.
    int           irecname {};     ///< Record name   1 : "ATOM  "; 2 : "HETATM".
    unsigned long pdbsegment {};   ///< Segment name.
    std::string   name;            ///< Pdb name.
    std::string   resname {""};    ///< Residue name.
    std::string   chainID {};      ///< Chain identifier.
    unsigned long resSeq {};       ///< Residue sequence number.
    float         occupancy {};    ///< Occupancy.
    float         tempFactor {};   ///< Temperature factor.
    std::string   element;         ///< Name of the element.
    std::string   charge;          ///< Electrostatic charge.
    float         vdWRad {};       ///< Van Der Waals radius.

    /// Constructor.
    explicit Pdb(
        const std::string& chainID,
        const std::string& resname,
        const std::string& atomname
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        uint ind,
        const std::string& chainID,
        uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3f& pos
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        uint ind,
        const std::string& chainID,
        uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3f& pos,
        float occupancy
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        uint ind,
        const std::string& chainID,
        uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3f& pos,
        float occupancy,
        float tempFactor
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        uint ind,
        const std::string& chainID,
        uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3f& pos,
        float occupancy,
        float tempFactor,
        const std::string& element
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        uint ind,
        const std::string& chainID,
        uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3f& pos,
        const std::string& element
    ) noexcept;

    /// Constructor.
    explicit Pdb(
        const std::string& record,
        unsigned long segm
    ) noexcept;

    /**
     * \brief Formated reading of a line record from a pdb file.
     * \param record Record in the pdb file.
     */
    void format_as_pdb(
        const std::string& record
    ) noexcept;

    /**
     * \brief Produces an atom-specific formatted record for a pdb file.
     * \param msgr Output message processor.
     */
    std::string format_as_pdb(
        Msgr& msgr
    ) const;

    /**
     * \brief Produces an atom-specific formatted record for a pdb file.
     * \param p Atom position coordinates.
     * \param msgr Output message processor.
     */
    std::string format_as_pdb(
        const A3f& p,
        Msgr& msgr
    ) const;

    /**
     * \brief Produces an atom-specific formatted record for a pdb file.
     * \param p Atom position coordinates.
     * \param msgr Output message processor.
     */
    std::string format_as_pdb(
        float* p,
        Msgr& msgr
    ) const;

    /**
     * \brief Produce an atom-specific formatted record for a pdb file.
     * \param a Atom-specific record.
     * \param msgr Output message processor.
     */
    static std::string format_as_pdb(
        const Pdb& a,
        Msgr& msgr
    );

    /// Produces an atom-specific formatted record for a pdb file.
    static std::string format_as_pdb(
        const std::string& recn,
        uint iatom,
        const std::string& atname,
        const std::string& resname,
        unsigned long ires,
        char ichain,
        const A3f& pos,
        float occ,
        float tempf,
        const std::string& elt,
        const std::string& chge,
        Msgr& msgr
    );

    /// Determines if the atom is part of an aminoacid.
    bool is_in_aminoacid() const noexcept;

    /**
     * \brief Reads from a pdb file.
     * \param file The file to read from.
     * \param a Atom-specific records.
     * \param other Non-atom records.
     * \param msgr Output message processor.
     */
    static void read(
        const std::filesystem::path& file,
        std::vector<Pdb>& a,
        std::vector<std::vector<std::string>>& other,
        Msgr& msgr
    );

    /**
     * \brief Writes to a pdb file.
     * \param file The file to write to.
     * \param a Atom-specific records.
     * \param other Non-atom records.
     * \param msgr Output message processor.
     */
    static void write(
        const std::filesystem::path& file,
        const std::vector<Pdb>& a,
        const std::vector<std::vector<std::string>>& other,
        Msgr& msgr
    );

private:

    /**
     * \brief Converts record name string to index.
     * \param s "HETATM" or "ATOM  ".
     */
    int irecn(const std::string& s) const noexcept;

    /**
     * \brief Gets record name.
     * \return "HETATM" or "ATOM  ".
     */
    std::string recname() const noexcept;

    /**
     * \brief Converts record name index to string.
     * \return "HETATM" or "ATOM  ".
     */
    static std::string irecname2str(int irecname) noexcept;

    /// Determines van Der Waals radius based on the atom data.
    float set_vdW() const noexcept;
};

}  // namespace utils::biochemical

#endif  // UTILS_BIOCHEMICAL_PDB_H
