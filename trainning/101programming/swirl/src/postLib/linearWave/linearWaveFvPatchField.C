/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "linearWaveFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    valueAbove_(pTraits<Type>::zero),
    valueBelow_(pTraits<Type>::zero),
    amplitude_(0),
    period_(1),
    curTimeIndex_(-1)
{}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    valueAbove_(pTraits<Type>(dict.lookup("valueAbove"))),
    valueBelow_(pTraits<Type>(dict.lookup("valueBelow"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    period_(readScalar(dict.lookup("period"))),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(valueBelow_);
    }
}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const linearWaveFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    valueAbove_(ptf.valueAbove_),
    valueBelow_(ptf.valueBelow_),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    curTimeIndex_(-1)
{}


template<class Type>
linearWaveFvPatchField<Type>::linearWaveFvPatchField
(
    const linearWaveFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    valueAbove_(ptf.valueAbove_),
    valueBelow_(ptf.valueBelow_),
    amplitude_(ptf.amplitude_),
    period_(ptf.period_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void linearWaveFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void linearWaveFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


// Update the coefficients associated with the patch field
template<class Type>
void linearWaveFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        scalar twoPi = 2* Foam::constant::mathematical::pi;

//         scalarField z = this->patch().Cf().component(vector::Z);
        scalarField z = this->patch().Cf().component(vector::Y);
        scalar h = amplitude_*
            Foam::sin(twoPi*this->db().time().value()/period_);

        Info << "t = " << this->db().time().value()
            << " h = " << h << endl;

        patchField = pos(z - h)*valueAbove_ + neg(z - h)*valueBelow_;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


// Write
template<class Type>
void linearWaveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("valueAbove")
        << valueAbove_ << token::END_STATEMENT << nl;
    os.writeKeyword("valueBelow")
        << valueBelow_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude")
        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("period")
        << period_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
