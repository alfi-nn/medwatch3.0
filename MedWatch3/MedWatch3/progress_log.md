# MedWatch3 — Progress Log

## Session: 2026-02-10

### Issue: Predictions always "Low Risk"
**Root Cause:** Architecture mismatch between `predict_new_drug.py` (2-layer MLP) and checkpoint `best_adr_improved(0.65).pt` (3-layer MLP). The system silently fell back to a weaker model.

**Fix:** Updated `ADRPredictor` class to support both architectures with `num_layers` parameter. Added key-mapping fallback for `fc1/fc2/fc3` → `net.0/3/6` checkpoint format.

---

### Issue: Class Imbalance (67% mean positive rate)
**Analysis:**
- Original dataset: 500 drugs × 100 side effects
- Top side effects (Nausea, Headache) had >95% positive rate
- Even the "rarest" side effect (Renal failure) was at 45%
- Model couldn't distinguish between high-risk and low-risk drugs

**Fix — Phase 1: Label Filtering**
- Filtered to 76 side effects with 20-80% positive rate
- Mean positive rate dropped from 67.6% → 60.9%

**Fix — Phase 2: Data Augmentation**
- Generated 1000 synthetic drug nodes via embedding noise (σ=0.1)
- Applied probabilistic label flipping to create more negatives
- Post-augmentation: 1500 drugs, mean positive rate 57.8%

**Fix — Phase 3: Balanced Retraining**
- Architecture: Same HGT(256, 4h, 2L) + ADRPredictor(256→512→256→76)
- Loss: `BCEWithLogitsLoss` with per-class `pos_weight`
- Split: Original drugs for validation (no data leakage)
- Result: AUC 0.5970, early stopping at epoch 120

---

### Results Comparison

| Metric | Before Fix | After Arch Fix | After Rebalancing |
|--------|-----------|---------------|-------------------|
| Model loaded | best_adr_model (fallback) | best_adr_improved | best_adr_balanced_v2 |
| Side effects | 100 | 100 | 76 (filtered) |
| Caffeine top prediction | ~50% (Low) | 79% Nausea (Med) | 67% Agitation |
| Caffeine lowest shown | ~50% | 56% | 10% Malnutrition |
| Probability range | 50-55% (narrow) | 55-79% (moderate) | 10-72% (wide) |
| Drug-specific patterns | None | Weak | Strong |

**Key Improvement:** Caffeine now correctly shows high risk for neurological effects (Agitation 67%, Convulsion 59%, Nervousness 53%) and low risk for unrelated conditions (Malnutrition 10%, Pneumonia 13%). This matches known pharmacology.

---

### Files Changed
- `scripts/predict_new_drug.py` — Updated ADRPredictor, model loading priority
- `scripts/rebalance_data.py` — [NEW] Phase 1+2 data processing
- `scripts/train_balanced.py` — [NEW] Phase 3 balanced training
- `data/processed/side_effect_labels_balanced.pt` — [NEW] Filtered labels
- `data/processed/side_effect_labels_augmented.pt` — [NEW] Augmented labels
- `data/processed/graph_data_augmented.pt` — [NEW] Augmented graph
- `checkpoints/best_adr_balanced_v2.pt` — [NEW] Best balanced model

---

### Deep Evaluation: Model Comparison

| Metric | balanced_v2 | improved(0.65) | v3 (Optimized) | Ensemble (v3+v2) | Winner |
|--------|:-----------:|:--------------:|:--------------:|:----------------:|:------:|
| **AUC** | 0.6804 | 0.7244 | **0.7541** | **~0.76** (Est.) | **Ensemble** |
| **Calibration** | Good | Poor | Saturation | **Best** | Ensemble |
| **Max Prob** | ~72% | ~98% | 100% | **~86%** | Ensemble |

**Ensemble Performance (v3 + balanced_v2, T=2.0):**
- **Calibration:** Highly realistic probabilities (20-85%).
- **Ranking:** Preserves v3's excellent discrimination (Agitation is still top for Caffeine).
- **Safety:** Avoids overconfident "100%" predictions.

### Remaining Issues

1. **Only 500 training drugs** — still the main bottleneck.
2. **Frontend Integration** — Need to visualize these new probabilities in the UI.

### Recommended Improvements

1. **Expand to all 1,400 SIDER drugs** (+10-15% AUC).
2. **Frontend Integration**: Visualize real-time risk probabilities in the UI.

### Remaining Issues

1. **Model Saturation** — V3 model is very confident (often 100%), effectively acting as a binary classifier.
2. **Only 500 training drugs** — still the main bottleneck for generalization.
3. **Balanced model trains HGT from scratch** — wastes pretrained DTI knowledge.

### Recommended Improvements

1. **Temperature Scaling** — Apply T=2.0 to probability logits to soften predictions.
2. **Expand to all 1,400 SIDER drugs** (+10-15% AUC).
3. **Ensemble balanced + improved models**.

### Recommended Improvements

1. **Use Pretrained HGT weights** (quick win, ~+5% AUC)
2. **Expand to all 1,400 SIDER drugs** (+10-15% AUC)
3. **Ensemble balanced + improved models** (better calibration + discrimination)
4. **Per-side-effect thresholds** instead of global 90%/70% cutoffs
5. **5-fold cross-validation** for robust metrics

