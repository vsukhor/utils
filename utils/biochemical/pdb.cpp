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

#include "pdb.h"

/// Library-wide.
namespace Utils {
/// \brief Biochemistry-related stuff.
namespace Biochemical {

float Pdb::scaling;    ///< Convert nm <-> A; nanometers are used internally.

Pdb::Pdb(
        const std::string& chainID,
        const std::string& resname,
        const std::string& atomname ) noexcept
    : name {atomname}
    , resname {resname}
    , chainID {chainID}
{}

Pdb::Pdb(
        const uint ind,
        const std::string& chainID,
        const uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3<float>& pos ) noexcept
    : pos {pos}
    , ind {ind}
    , name {atomname}
    , resname {resname}
    , chainID {chainID}
    , resSeq {resSeq}
    , occupancy {}
    , tempFactor {}
    , element {}
    , charge {}
{}

Pdb::Pdb(
        const uint ind,
        const std::string& chainID,
        const uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3<float>& pos,
        const float occupancy ) noexcept
    : pos {pos}
    , ind {ind}
    , name {atomname}
    , resname {resname}
    , chainID {chainID}
    , resSeq {resSeq}
    , occupancy {occupancy}
{}

Pdb::Pdb(
        const uint ind,
        const std::string& chainID,
        const uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3<float>& pos,
        const float occupancy,
        const float tempFactor ) noexcept
    : pos {pos}
    , ind {ind}
    , name {atomname}
    , resname {resname}
    , chainID {chainID}
    , resSeq {resSeq}
    , occupancy {occupancy}
    , tempFactor {tempFactor}
{}

Pdb::Pdb(
        const uint ind,
        const std::string& chainID,
        const uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3<float>& pos,
        const float occupancy,
        const float tempFactor,
        const std::string& element ) noexcept
    : pos {pos}
    , ind {ind}
    , name {atomname}
    , resname {resname}
    , chainID {chainID}
    , resSeq {resSeq}
    , occupancy {occupancy}
    , tempFactor {tempFactor}
    , element {element}
{}

Pdb::Pdb(
        const uint ind,
        const std::string& chainID,
        const uint resSeq,
        const std::string& resname,
        const std::string& atomname,
        const A3<float>& pos,
        const std::string& element ) noexcept
    : pos {pos}
    , ind {ind}
    , name {atomname}
    , resname {resname}
    , chainID {chainID}
    , resSeq {resSeq}
    , element {element}
{}

Pdb::Pdb(
        const std::string record,
        const ulong segm ) noexcept
    : pdbsegment(segm) 
{
    format_as_pdb(record);
}

void Pdb::
format_as_pdb( const std::string& record ) noexcept
{
    irecname = irecn(record.substr(0, 6));
    ind = static_cast<uint>(std::stoi(trim(record.substr(6, 5))));
    name = record.substr(12, 4);
    resname = record.substr(17, 3);
    chainID = record[21];
    resSeq = std::stoi(trim(record.substr(22, 4)));

    // Converts to nm from data read in A:
    pos[0] = scaling * static_cast<float>(std::stod(trim(record.substr(30, 8))));
    pos[1] = scaling * static_cast<float>(std::stod(trim(record.substr(38, 8))));
    pos[2] = scaling * static_cast<float>(std::stod(trim(record.substr(46, 8))));

    occupancy = static_cast<float>(std::stod(trim(record.substr(54, 6))));
    tempFactor = static_cast<float>(std::stod(trim(record.substr(60, 6))));
    element = record.substr(76, 2);
    charge  = record.substr(78, 2);
    vdWRad = set_vdW();
}

std::string Pdb::
format_as_pdb( Msgr& msgr ) const
{
    return format_as_pdb(recname(), ind, name, resname, resSeq, chainID[0],
                         pos, occupancy, tempFactor, element, charge, msgr);
}

std::string Pdb::
format_as_pdb( const A3<float>& p, Msgr& msgr ) const
{
    return format_as_pdb(recname(), ind, name, resname, resSeq, chainID[0],
                         p, occupancy, tempFactor, element, charge, msgr);
}

std::string Pdb::
format_as_pdb( float* p, Msgr& msgr ) const
{
    return format_as_pdb(recname(), ind, name, resname, resSeq, chainID[0],
                         {p[0], p[1], p[2]}, occupancy, tempFactor, element,
                         charge, msgr);
}

std::string Pdb::
format_as_pdb( const Pdb& a, Msgr& msgr )
{
    return format_as_pdb(irecname2str(a.irecname), a.ind, a.name, a.resname,
                         a.resSeq, a.chainID.c_str()[0], a.pos, a.occupancy,
                         a.tempFactor, a.element, a.charge, msgr);
}

std::string Pdb::
format_as_pdb(
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
    Msgr& msgr )
{
    // http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    std::ostringstream record;

    // columns 1-6  Record name  "ATOM  " or "TER   "
    record << recn;
    record.width(5);

    record << std::right;

    // columns 7-11   Integer   Pdb serial number
    record << STR(iatom);

    // column  12  empty
    record << " ";

    // columns 13-16  Pdb name  Pdb name
    record.width(4);
    record << atname;

    // column 17   Character   altLoc   Alternate location indicator
    record << " ";

    // columns 18-20  Residue name  resName   Residue name
    record.width(3);
    record << resname;

    // column 21  empty
    record << " ";

    // column 22   Character  chainID Chain identifier
    record.width(1);
    record << ichain;

    // columns 23-26  Integer  resSeq   Residue sequence number
    record.width(4);
    record << STR(uint(ires));

    // column  27   AChar  iCode  Code for insertion of residues
    record << " ";

    if (recn == "ATOM  " || recn == "HETATM") {
        A3<float> p {pos/scaling};
        float mx{1000.f};
        if (std::abs(p[0]) >= mx ||
            std::abs(p[1]) >= mx ||
            std::abs(p[2]) >= mx) {
            msgr.exit("Error: Pdb position exceeds pdb numerical limit ");
        }

        // columns 28-30      empty
        record << "   ";

        record.setf(std::ios::fixed);

        // columns 31-38  Real(8.3)  x  Orthogonal coordinates for X in Angstroms
        record.width(8);
        record.precision(3);
        record << p[0];

        // columns 39-46  Real(8.3)  y  Orthogonal coordinates for Y in Angstroms
        record.width(8);
        record.precision(3);
        record << p[1];

        // columns 47-54  Real(8.3)  z  Orthogonal coordinates for Z in Angstroms
        record.width(8);
        record.precision(3);
        record << p[2];

        // columns 55-60   Real(6.2)  occupancy   Occupancy
        record.width(6);
        record.precision(2);
        record << occ;

        // columns 61-66  Real(6.2)  tempFactor  Temperature  factor
        record.width(6);
        record.precision(2);
        record << tempf;

        // columns 67-76    empty
        record << "          ";

        // columns 77-78  LString(2)  element  Element symbol, right-justified
        record << std::right;
        record.width(2);
        record << elt;

        // columns 79-80  LString(2)  charge  Charge  on the atom (e.g. "2+" or "1-")
        record.width(2);
        record << chge;
    }

    return record.str();
}

int Pdb::
irecn( const std::string& s ) const noexcept
{
    return (s == "ATOM  ") ? 1
                           : ((s == "HETATM") ? 2 : -1);
} 

std::string Pdb::
recname() const noexcept
{
    return (irecname == 1) ? "ATOM  "
                           : "HETATM";
} 

std::string Pdb::
irecname2str( const int irecname ) noexcept
{
    return (irecname == 1) ? "ATOM  "
                           : "HETATM";
}

float Pdb::
set_vdW() const noexcept
{
    if (resname == "UNX" || resname == "UNK")
        return 0.4f;
    else
        switch (name.front()) {
            case 'C': return 0.170f;
            case 'N': return 0.155f;
            case 'O': return 0.152f;
            case 'H': return 0.120f;
            case 'S': return 0.180f;
            default:  return 0.160f;
        }
}

bool Pdb::
is_in_aminoacid() const noexcept
{
    std::vector<std::string> resnames { "ILE", "LEU", "PHE", "VAL", "MET", "PRO", "TRP", "HIS", "THR", "GLU",
                                        "GLN", "CYS", "TYR", "ALA", "SER", "ASN", "ASP", "ARG", "GLY", "LYS" };
    return std::find(resnames.begin(), resnames.end(), resname) != resnames.end();
}

void Pdb::
read( const std::string& filename,
      std::vector<Pdb>& a,
      vec2<std::string>& other,
      Msgr& msgr )
{
    std::ifstream fin {filename};
    if (!fin.is_open())
        msgr.exit("Unable to open file for reading at "+filename);

    a.clear();
    other.clear();
    std::string line;
    bool newr {true};
    while (std::getline(fin, line)) {
        const auto s {line.substr(0, 6)};
        if (s == "ATOM  " || s == "HETATM") {
            a.emplace_back(line, other.size());
            newr = true;
        }
        else {
            if (newr) {
                other.push_back(std::vector<std::string>());
                newr = false;
            }
            other.back().push_back(line);
        }
    }
}

void Pdb::
write( const std::string& filename,
       const std::vector<Pdb>& a,
       const vec2<std::string>& other,
       Msgr& msgr )
{
    std::ofstream fout {filename};
    if (!fout.is_open())
        msgr.exit("Unable to open file for writing at "+filename);

    szt i = 0;
    if (!other.empty())
        for (szt j=0; j<other.size(); j++) {
            while (i < a.size() && a[i].pdbsegment==j) {
                fout << a[i].format_as_pdb(msgr) << std::endl;
                i++;
            }
            for (szt k=0; k<other[j].size(); k++)
                fout << other[j][k] << std::endl;
        }
    else
        while (i < a.size()) {
            fout << a[i].format_as_pdb(msgr) << std::endl;
            i++;
        }
}

}    // namespace Biochemical
}    // namespace Utils








