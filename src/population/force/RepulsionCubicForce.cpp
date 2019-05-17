#include "RepulsionCubicForce.hpp"

template<unsigned DIM>
RepulsionCubicForce<DIM>::RepulsionCubicForce()
   : CubicGeneralisedLinearSpringForce<DIM>()
{
}

template<unsigned DIM>
void RepulsionCubicForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a NodeBasedCellPopulation
    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("RepulsionCubicForce is to be used with a NodeBasedCellPopulation only");
    }

    std::vector< std::pair<Node<DIM>*, Node<DIM>* > >& r_node_pairs = (static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetNodePairs();

    for (typename std::vector< std::pair<Node<DIM>*, Node<DIM>* > >::iterator iter = r_node_pairs.begin();
        iter != r_node_pairs.end();
        iter++)
    {
        std::pair<Node<DIM>*, Node<DIM>* > pair = *iter;

        Node<DIM>* p_node_a = pair.first;
        Node<DIM>* p_node_b = pair.second;

        // Get the node locations
        const c_vector<double, DIM>& r_node_a_location = p_node_a->rGetLocation();
        const c_vector<double, DIM>& r_node_b_location = p_node_b->rGetLocation();

        // Get the node radii
        double node_a_radius = p_node_a->GetRadius();
        double node_b_radius = p_node_b->GetRadius();

        // Get the unit vector parallel to the line joining the two nodes
        c_vector<double, DIM> unit_difference;

        unit_difference = (static_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))->rGetMesh().GetVectorFromAtoB(r_node_a_location, r_node_b_location);

        // Calculate the value of the rest length
        double rest_length = node_a_radius+node_b_radius;

        if (norm_2(unit_difference) < rest_length)
        {
            // Calculate the force between nodes
            c_vector<double, DIM> force = this->CalculateForceBetweenNodes(p_node_a->GetIndex(), p_node_b->GetIndex(), rCellPopulation);
            c_vector<double, DIM> negative_force = -1.0 * force;
            for (unsigned j=0; j<DIM; j++)
            {
                assert(!std::isnan(force[j]));
            }
            // Add the force contribution to each node
            p_node_a->AddAppliedForceContribution(force);
            p_node_b->AddAppliedForceContribution(negative_force);
        }
    }
}

template<unsigned DIM>
void RepulsionCubicForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    CubicGeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RepulsionCubicForce<1>;
template class RepulsionCubicForce<2>;
template class RepulsionCubicForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionCubicForce)
