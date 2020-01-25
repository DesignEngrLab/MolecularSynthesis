using System.Linq;
using System.Xml.Schema;
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
            var result = new double[10];
            // call openBabel functions to get x, y, z's into
            
            foreach (var n in cand.graph.nodes)
            {
                var mass = 0.0;// OBFunctions.
                result[0] += mass;
                result[1] += mass * n.X;
                result[2] += mass * n.Y;
                result[3] += mass * n.Z;
                result[1] += mass * n.X;
                result[1] += mass * n.X;
                result[1] += mass * n.X;
            }

            return result;
        }




    }
}