/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gaussGrad.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
gaussGrad<Type>::gradf
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        new GeometricField<GradType, fvPatchField, volMesh>
        (
            IOobject
            (
                name,
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<GradType>
            (
                "0",
                ssf.dimensions()/dimLength,
                pTraits<GradType>::zero
            ),
            zeroGradientFvPatchField<GradType>::typeName
        )
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf();



    Field<GradType>& igGrad = gGrad;
    const Field<Type>& issf = ssf;

    //Info << "ssf = " << ssf << endl;

    forAll(owner, facei)
    {
        GradType Sfssf = Sf[facei]*issf[facei];

        igGrad[owner[facei]] += Sfssf;
        igGrad[neighbour[facei]] -= Sfssf;

        //Info << facei << "igGrad = " << igGrad << endl;
    }

    //Info << "igGrad = " << igGrad << endl;

    forAll(mesh.boundary(), patchi)
    {
        const unallocLabelList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        const fvsPatchField<Type>& pssf = ssf.boundaryField()[patchi];
        //Info << "patch type = " << mesh.boundary()[patchi].patch().type() << endl;

        if (mesh.boundary()[patchi].patch().type() != "processor")
        {
            //Info << "Went inside" << endl;
            forAll(mesh.boundary()[patchi], facei)
            {
                //Info << "pFaceCells = " << pFaceCells[facei] << endl;
                igGrad[pFaceCells[facei]] += pSf[facei]*pssf[facei];    
            
            }
        //Info << "igGrad = " << igGrad << endl;

        }
        
        //if (mesh.boundaryField()[patchi].patch().type() == "processor")
        
        //Info << "igGrad = " << igGrad << endl;
    }

    //Info << "V = " << mesh.V() << endl;
    igGrad /= mesh.V();

    //Info << "igGrad = " << igGrad << endl;

    //Info << "in gradf " << endl;

    gGrad.correctBoundaryConditions();
    
    //Info << "igGrad after correct = " << igGrad << endl;

    return tgGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
gaussGrad<Type>::calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    //Info << "vsf = " << vsf << endl;

    //Info << "interpolate = " << tinterpScheme_().interpolate(vsf) << endl;

    tmp<GeometricField<GradType, fvPatchField, volMesh> > tgGrad
    (
        gradf(tinterpScheme_().interpolate(vsf), name)
    );
    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad();

    gGrad.rename("grad(" + vsf.name() + ')');
    //Info << "calc grad correct BC " << endl;
    this->correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}

template<class Type>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
> gaussGrad<Type>::fvmGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> fvmGrad\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit gradient operator defined only for scalar."
        << abort(FatalError);

    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<BlockLduSystem<vector, GradType> > tbs
    (
        new BlockLduSystem<vector, GradType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
