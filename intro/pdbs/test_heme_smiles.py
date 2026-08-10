from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
candidates={
'cand1':'Cc1c2[n-]3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe+2]36[N]7=C(C=C8[N-]6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O',
'cand2':'Cc1c2[n-]3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe+2]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O',
'neutral_rcsb':'Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O'
}
for name,s in candidates.items():
 m=Chem.MolFromSmiles(s)
 print(name,'parse',m is not None)
 if m:
  print(' canonical',Chem.MolToSmiles(m,canonical=True))
  print(' charge',Chem.GetFormalCharge(m),'atoms',m.GetNumAtoms(),'heavy',m.GetNumHeavyAtoms(),'formula',rdMolDescriptors.CalcMolFormula(m))
  print(' metal',[(a.GetIdx(),a.GetSymbol(),a.GetFormalCharge(),a.GetDegree()) for a in m.GetAtoms() if a.GetAtomicNum()==26])
  print(' chargedN',[(a.GetIdx(),a.GetSymbol(),a.GetFormalCharge(),a.GetIsAromatic()) for a in m.GetAtoms() if a.GetAtomicNum()==7 and a.GetFormalCharge()!=0])
