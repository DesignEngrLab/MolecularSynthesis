using GraphSynth.Representation;
using System;
using System.Collections.Generic;
using System.Data;
using System.Text;

namespace MolecularSynthesis.GS.Plugin
{
    public class TreeCandidate : candidate
    {
        public TreeCandidate Parent; //as a field
        public List<TreeCandidate> Children;
        public double S;
        public double n;
        public double UCB;

        public TreeCandidate(candidate seedCandidate) 
        {
            this.activeRuleSetIndex = seedCandidate.activeRuleSetIndex;
            this.age = seedCandidate.age;
            this.designParameters = seedCandidate.designParameters;
            this.GenerationStatus = seedCandidate.GenerationStatus;
            this.graph = seedCandidate.graph;
            this.graphFileName = seedCandidate.graphFileName;
            this.performanceParams = seedCandidate.performanceParams;
        }
        public override candidate copy()
        {
            return new TreeCandidate(base.copy());
        }
    }
}
