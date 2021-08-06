

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
            this.recipe = seedCandidate.recipe;
        }
        public override candidate copy()
        {
            var copytreecandidate= new TreeCandidate(base.copy());
            copytreecandidate.Parent = this.Parent;
            copytreecandidate.Children = new List<TreeCandidate>();
            copytreecandidate.S = S;
            copytreecandidate.n = n;
            copytreecandidate.UCB = UCB;

            //copytreecandidate.recipe = this.recipe;
            //copytreecandidate.graph = this.graph;




            return copytreecandidate;
        }

        public TreeCandidate deepcopy()
        {
            TreeCandidate other = (TreeCandidate)this.MemberwiseClone();
            //other.Parent = new TreeCandidate(IdInfo.IdNumber);
            
            return other;

            
        }
             
    }
}
