"""Frozen additive ipSAE2/water consumer, independent of live upstream publishers."""
import json
import pandas as pd

def configure(store):
 data=store.data
 store.extra_metrics=store.read('enrichment_metrics.json',[])
 store.metric_labels=dict(store.metric_labels,**{r['id']:r['label'] for r in store.extra_metrics})
 store.extra_scores=pd.read_parquet(data/'enrichment_scores.parquet').set_index('pose_id') if (data/'enrichment_scores.parquet').exists() else pd.DataFrame()
 if len(store.extra_scores):
  # Outer additive join; wet and dry native IDs cannot collide.
  store.scores=store.scores.join(store.extra_scores,how='outer')
 store.water_energy=pd.read_parquet(data/'enrichment_water_physics.parquet').set_index('pose_id').to_dict('index') if (data/'enrichment_water_physics.parquet').exists() else {}
 for r in store.read('enrichment_ligands.json',[]):store.ligands[r['id']]=r
 for r in store.read('enrichment_ledger.json',[]):
  assert r['pose_id'] not in store.ledger
  store.ledger[r['pose_id']]=r
 store.receptor.update(store.read('enrichment_receptor.json',{}))
 store.enrichment=store.read('enrichment_build.json',{})

def decorate(store,pose,pid):
 if pid not in store.scores.index:return
 row=store.scores.loc[pid]
 pose['native_confidence_percentile']=row.get('native_confidence_percentile')
 pose['native_confidence_missing_reason']=row.get('native_confidence_missing_reason')
 pose['native_confidence_score']=row.get('native_confidence_score')
 pose['native_confidence_score_kind']=row.get('native_confidence_score_kind')
 pose['native_confidence_pair_semantics']=row.get('native_confidence_pair_semantics')
 pose['native_confidence_reference_id']=row.get('native_confidence_reference_id')
 pose['arm']=row.get('arm','dry')
 pose['source_pose_id']=row.get('source_pose_id',pid)
 if pose['arm']=='water':pose['scores']['pl_iptm']=row.get('native_confidence_score')
 # Raw encodings are source values, not locally fitted percentile substitutes.
 for m in store.extra_metrics:
  # Only the actual percentile is a percentile. ipSAE2 remains raw numeric.
  if m['id']=='native_confidence_percentile':pose['percentiles'][m['id']]=row.get(m['id'])
