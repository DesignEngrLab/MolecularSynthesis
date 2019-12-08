using GraphSynth.Representation;
using StarMathLib;



namespace PropertyEvaluation
{
    public class Evaluation
    {
        public static int CountAtoms(candidate cand)
        {
            return cand.graph.nodes.Count;
        }

        public static double[] CalcMoment(candidate cand)
        {
            return new double[]{.0,.0,.0};
        }




    }
}