#include "PatchSizeTrackingModifier.hpp"

template<unsigned DIM>
PatchSizeTrackingModifier<DIM>::PatchSizeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PatchSizeTrackingModifier<DIM>::~PatchSizeTrackingModifier()
{
}

template<unsigned DIM>
void PatchSizeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PatchSizeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PatchSizeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    std::map<int, unsigned> clonalPatches;

    // First calculate the patch size for each ancestor
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get ancestor ID (set to -1 if unset)
        int ancestor = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1) : cell_iter->GetAncestor();
        if (clonalPatches.count(ancestor)){
            ++clonalPatches[ancestor]; // increase the count
        }
        else{
            clonalPatches[ancestor] = 1 ;
        }
                
    }

    // Next iterate over the population to store each cell's patch size
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get ancestor ID (set to -1 if unset)
        int ancestor = (cell_iter->GetAncestor() == UNSIGNED_UNSET) ? (-1) : cell_iter->GetAncestor();
        cell_iter->GetCellData()->SetItem("patch size", clonalPatches[ancestor]);
    }
}

template<unsigned DIM>
void PatchSizeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PatchSizeTrackingModifier<1>;
template class PatchSizeTrackingModifier<2>;
template class PatchSizeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PatchSizeTrackingModifier)