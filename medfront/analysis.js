/* ================================================================
   ADR-Predict Dashboard — JavaScript Logic
   ================================================================ */

const API_BASE = 'http://localhost:8000';
let history = [];
let historyIdx = -1;

// Organ mapping for radar chart
const ORGAN_MAP = {
    'Liver': ['Hepatitis', 'Hepatotoxicity', 'Hepatobiliary disease'],
    'Heart': ['Tachycardia', 'Arrhythmia', 'Cardiac disorder', 'Myocardial infarction', 'Palpitations'],
    'Kidney': ['Renal failure', 'Haematuria', 'Dysuria', 'Urinary tract disorder', 'Pollakiuria', 'Urinary tract infection'],
    'CNS': ['Tremor', 'Convulsion', 'Confusional state', 'Somnolence', 'Nervousness', 'Agitation', 'Anxiety', 'Mental disorder', 'Nervous system disorder', 'Hypoaesthesia', 'Neuropathy peripheral'],
    'GI': ['Flatulence', 'Stomatitis', 'Gastrointestinal disorder', 'Dysgeusia'],
    'Skin': ['Erythema', 'Alopecia', 'Erythema multiforme', 'Angioedema', 'Skin disorder', 'Hyperhidrosis', 'Flushing'],
};

// Medical icons per side effect
const SE_ICONS = {
    'Hepatitis': '🫀', 'Hepatotoxicity': '🫀', 'Hepatobiliary disease': '🫀',
    'Tachycardia': '❤️', 'Arrhythmia': '❤️', 'Cardiac disorder': '❤️', 'Myocardial infarction': '❤️', 'Palpitations': '❤️',
    'Renal failure': '🫘', 'Haematuria': '🫘', 'Dysuria': '🫘',
    'Tremor': '🧠', 'Convulsion': '🧠', 'Confusional state': '🧠', 'Somnolence': '🧠',
    'Alopecia': '💇', 'Erythema': '🔬', 'Nausea': '🤢',
};

/* ---------- DOM ---------- */
const $ = (s) => document.querySelector(s);
const $$ = (s) => document.querySelectorAll(s);

const smilesInput = $('#smiles-input');
const analyzeBtn = $('#analyze-btn');
const clearBtn = $('#clear-btn');
const pasteBtn = $('#paste-btn');

const loadingState = $('#loading-state');
const emptyState = $('#empty-state');
const seCards = $('#se-cards');
const chartsSection = $('#charts-section');
const molPlaceholder = $('#mol-placeholder');
const molActive = $('#mol-active');
const riskOverlay = $('#risk-overlay');

/* ---------- INIT ---------- */
document.addEventListener('DOMContentLoaded', () => {
    analyzeBtn.addEventListener('click', handleAnalyze);
    clearBtn.addEventListener('click', handleClear);
    pasteBtn.addEventListener('click', handlePaste);
    smilesInput.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') handleAnalyze();
    });

    // Filter tabs
    $$('.filter-tab').forEach(tab => {
        tab.addEventListener('click', () => {
            $$('.filter-tab').forEach(t => t.classList.remove('active'));
            tab.classList.add('active');
        });
    });

    // Nav icons
    $$('.topnav__icon-btn').forEach(btn => {
        btn.addEventListener('click', () => {
            $$('.topnav__icon-btn').forEach(b => b.classList.remove('active'));
            btn.classList.add('active');
        });
    });

    // Quick-start buttons
    $$('.quick-start__btn').forEach(btn => {
        btn.addEventListener('click', () => {
            smilesInput.value = btn.dataset.smiles;
            handleAnalyze();
        });
    });

    // History nav
    $('#history-prev').addEventListener('click', () => navigateHistory(-1));
    $('#history-next').addEventListener('click', () => navigateHistory(1));
});

/* ---------- HANDLE ANALYZE ---------- */
async function handleAnalyze() {
    const smiles = smilesInput.value.trim();
    if (!smiles) {
        smilesInput.focus();
        smilesInput.style.borderColor = '#EF4444';
        setTimeout(() => smilesInput.style.borderColor = '', 1500);
        return;
    }

    // Show loading
    emptyState.classList.add('hidden');
    seCards.classList.add('hidden');
    loadingState.classList.remove('hidden');
    chartsSection.classList.add('hidden');
    riskOverlay.classList.add('hidden');
    analyzeBtn.classList.add('loading');
    analyzeBtn.innerHTML = `<span class="loading-spinner" style="width:16px;height:16px;border-width:2px;"></span> Analyzing…`;

    const startTime = performance.now();

    try {
        const resp = await fetch(`${API_BASE}/predict`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles }),
        });

        if (!resp.ok) {
            const err = await resp.json();
            throw new Error(err.detail || `HTTP ${resp.status}`);
        }

        const data = await resp.json();
        const latency = ((performance.now() - startTime) / 1000).toFixed(1);

        // Update status
        $('#confidence-val').textContent = 'High';
        $('#latency-val').textContent = `${latency}s`;
        $('#model-val').textContent = 'Ensemble GNN + Transformer';

        // Process results
        const results = data.results;
        renderPredictions(results);
        renderBarChart(results);
        renderRadarChart(results);
        renderMolecule(smiles, results);
        updateRiskOverlay(results);
        updateFooter(smiles);

        // Add to history
        history.push({ smiles, results, timestamp: new Date() });
        historyIdx = history.length - 1;
        updateHistoryNav();

    } catch (err) {
        loadingState.classList.add('hidden');
        emptyState.classList.remove('hidden');
        emptyState.querySelector('p').innerHTML = `<span style="color:var(--red)">Error: ${err.message}</span><br/>Check that the API server is running on port 8000.`;
    } finally {
        analyzeBtn.classList.remove('loading');
        analyzeBtn.innerHTML = `<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="22 12 18 12 15 21 9 3 6 12 2 12"/></svg> Analyze Molecule`;
    }
}

/* ---------- RENDER PREDICTIONS ---------- */
function renderPredictions(results) {
    loadingState.classList.add('hidden');
    seCards.classList.remove('hidden');
    seCards.innerHTML = '';

    // Sort by probability descending
    const sorted = Object.entries(results).sort((a, b) => b[1] - a[1]);

    // Show top 20 for cleanliness
    const top = sorted.slice(0, 20);

    top.forEach(([name, prob], i) => {
        const pct = (prob * 100).toFixed(1);
        const severity = prob > 0.7 ? 'high' : prob > 0.4 ? 'moderate' : 'low';
        const icon = SE_ICONS[name] || '💊';
        const barColor = getGradientColor(prob);

        const card = document.createElement('div');
        card.className = 'se-card';
        card.style.animationDelay = `${i * 0.04}s`;
        card.innerHTML = `
      <div class="se-card__icon se-card__icon--${severity}">${icon}</div>
      <div class="se-card__info">
        <div class="se-card__name">${name}</div>
        <div class="se-card__bar-wrap">
          <div class="se-card__bar" style="width:0%;background:${barColor}"></div>
        </div>
      </div>
      <div class="se-card__meta">
        <span class="se-card__pct se-card__pct--${severity}">${pct}%</span>
        <span class="se-card__severity se-card__severity--${severity}">${severity}</span>
      </div>
    `;
        seCards.appendChild(card);

        // Animate bars
        requestAnimationFrame(() => {
            setTimeout(() => {
                card.querySelector('.se-card__bar').style.width = `${pct}%`;
            }, 100 + i * 40);
        });
    });
}

/* ---------- BAR CHART ---------- */
function renderBarChart(results) {
    chartsSection.classList.remove('hidden');
    const chart = $('#bar-chart');
    chart.innerHTML = '';

    const sorted = Object.entries(results).sort((a, b) => b[1] - a[1]);
    const top10 = sorted.slice(0, 10);

    top10.forEach(([name, prob], i) => {
        const pct = (prob * 100).toFixed(1);
        const color = getGradientColor(prob);
        const row = document.createElement('div');
        row.className = 'bar-row';
        const isSmall = prob < 0.12;
        row.innerHTML = `
      <span class="bar-row__label">${name}</span>
      <div class="bar-row__track">
        <div class="bar-row__fill" style="width:0%;background:${color}">
          ${!isSmall ? `<span class="bar-row__pct">${pct}%</span>` : ''}
        </div>
      </div>
      ${isSmall ? `<span class="bar-row__pct--outside">${pct}%</span>` : ''}
    `;
        chart.appendChild(row);

        requestAnimationFrame(() => {
            setTimeout(() => {
                row.querySelector('.bar-row__fill').style.width = `${pct}%`;
            }, 200 + i * 80);
        });
    });
}

/* ---------- RADAR CHART ---------- */
function renderRadarChart(results) {
    const canvas = $('#radar-canvas');
    const ctx = canvas.getContext('2d');
    const W = canvas.width;
    const H = canvas.height;
    const cx = W / 2;
    const cy = H / 2;
    const maxR = Math.min(cx, cy) - 50;

    ctx.clearRect(0, 0, W, H);

    const organs = Object.keys(ORGAN_MAP);
    const n = organs.length;
    const angleStep = (Math.PI * 2) / n;

    // Compute organ scores
    const scores = organs.map(organ => {
        const keys = ORGAN_MAP[organ];
        let sum = 0, count = 0;
        for (const [name, prob] of Object.entries(results)) {
            if (keys.some(k => name.toLowerCase().includes(k.toLowerCase()))) {
                sum += prob; count++;
            }
        }
        return count > 0 ? sum / count : 0;
    });

    // Draw grid
    [0.2, 0.4, 0.6, 0.8, 1.0].forEach(level => {
        ctx.beginPath();
        for (let i = 0; i <= n; i++) {
            const angle = i * angleStep - Math.PI / 2;
            const x = cx + Math.cos(angle) * maxR * level;
            const y = cy + Math.sin(angle) * maxR * level;
            if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.closePath();
        ctx.strokeStyle = '#E8ECF0';
        ctx.lineWidth = 1;
        ctx.stroke();
    });

    // Draw axes
    for (let i = 0; i < n; i++) {
        const angle = i * angleStep - Math.PI / 2;
        ctx.beginPath();
        ctx.moveTo(cx, cy);
        ctx.lineTo(cx + Math.cos(angle) * maxR, cy + Math.sin(angle) * maxR);
        ctx.strokeStyle = '#E8ECF0';
        ctx.stroke();
    }

    // Draw labels
    ctx.font = '500 11px Inter, sans-serif';
    ctx.fillStyle = '#6B7280';
    ctx.textAlign = 'center';
    for (let i = 0; i < n; i++) {
        const angle = i * angleStep - Math.PI / 2;
        const lx = cx + Math.cos(angle) * (maxR + 28);
        const ly = cy + Math.sin(angle) * (maxR + 28);
        ctx.fillText(organs[i], lx, ly + 4);
    }

    // Draw data polygon
    ctx.beginPath();
    for (let i = 0; i <= n; i++) {
        const idx = i % n;
        const angle = idx * angleStep - Math.PI / 2;
        const r = scores[idx] * maxR;
        const x = cx + Math.cos(angle) * r;
        const y = cy + Math.sin(angle) * r;
        if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.closePath();
    ctx.fillStyle = 'rgba(13, 148, 136, 0.12)';
    ctx.fill();
    ctx.strokeStyle = '#0D9488';
    ctx.lineWidth = 2;
    ctx.stroke();

    // Draw data points
    for (let i = 0; i < n; i++) {
        const angle = i * angleStep - Math.PI / 2;
        const r = scores[i] * maxR;
        const x = cx + Math.cos(angle) * r;
        const y = cy + Math.sin(angle) * r;
        ctx.beginPath();
        ctx.arc(x, y, 4, 0, Math.PI * 2);
        ctx.fillStyle = '#0D9488';
        ctx.fill();
        ctx.strokeStyle = 'white';
        ctx.lineWidth = 2;
        ctx.stroke();
    }
}

/* ---------- MOLECULE VIZ ---------- */
function renderMolecule(smiles, results) {
    molPlaceholder.classList.add('hidden');
    molActive.classList.remove('hidden');
    // Risk info is now embedded in the result card, no need for the overlay
    riskOverlay.classList.add('hidden');

    const container = molActive;
    container.innerHTML = '';

    // Calculate summary stats from results
    const probs = Object.values(results || {});
    const avg = probs.length > 0 ? probs.reduce((a, b) => a + b, 0) / probs.length : 0;
    const highCount = probs.filter(p => p > 0.6).length;
    const modCount = probs.filter(p => p > 0.35 && p <= 0.6).length;
    const lowCount = probs.filter(p => p <= 0.35).length;
    const maxProb = probs.length > 0 ? Math.max(...probs) : 0;
    const riskPct = (avg * 100).toFixed(1);

    // Risk level
    const riskLevel = avg > 0.6 ? 'High' : avg > 0.35 ? 'Moderate' : 'Low';
    const riskColor = avg > 0.6 ? '#EF4444' : avg > 0.35 ? '#F59E0B' : '#10B981';
    const riskColorClass = avg > 0.6 ? 'high' : avg > 0.35 ? 'moderate' : 'low';

    // Build gauge SVG — semi-circle arc
    const gaugeAngle = avg * 180; // 0-180 degrees
    const gaugeEndX = 100 - 80 * Math.cos(gaugeAngle * Math.PI / 180);
    const gaugeEndY = 100 - 80 * Math.sin(gaugeAngle * Math.PI / 180);
    const largeArc = gaugeAngle > 90 ? 1 : 0;

    const card = document.createElement('div');
    card.className = 'mol-result-card';
    card.innerHTML = `
        <div class="mol-result-card__title">Analysis Complete</div>
        <div class="mol-result-card__subtitle">Adverse Drug Reaction Prediction Summary</div>
        <div class="mol-result-card__smiles">${smiles}</div>

        <div class="risk-gauge">
          <svg class="risk-gauge__svg" viewBox="0 0 200 110">
            <!-- Background arc -->
            <path d="M 20 100 A 80 80 0 0 1 180 100" fill="none" stroke="#E8ECF0" stroke-width="12" stroke-linecap="round"/>
            <!-- Colored arc -->
            <path d="M 20 100 A 80 80 0 ${largeArc} 1 ${gaugeEndX.toFixed(1)} ${gaugeEndY.toFixed(1)}" fill="none" stroke="${riskColor}" stroke-width="12" stroke-linecap="round" style="transition: all 1s ease"/>
            <!-- Center text -->
            <text x="100" y="88" text-anchor="middle" fill="${riskColor}" font-size="28" font-weight="800" font-family="Inter, sans-serif">${riskPct}%</text>
            <text x="100" y="106" text-anchor="middle" fill="#9CA3AF" font-size="10" font-weight="600" font-family="Inter, sans-serif">${riskLevel} Risk</text>
          </svg>
        </div>

        <div class="mol-result-stats">
          <div class="mol-result-stat">
            <span class="mol-result-stat__value mol-result-stat__value--high">${highCount}</span>
            <span class="mol-result-stat__label">High Risk</span>
          </div>
          <div class="mol-result-stat">
            <span class="mol-result-stat__value mol-result-stat__value--moderate">${modCount}</span>
            <span class="mol-result-stat__label">Moderate</span>
          </div>
          <div class="mol-result-stat">
            <span class="mol-result-stat__value mol-result-stat__value--low">${lowCount}</span>
            <span class="mol-result-stat__label">Low Risk</span>
          </div>
        </div>

        <div style="font-size:0.72rem; color:#9CA3AF; margin-top:4px;">
          ${probs.length} side effects analyzed · Peak probability: ${(maxProb * 100).toFixed(1)}%
        </div>
    `;

    container.appendChild(card);
}

/* ---------- RISK OVERLAY ---------- */
function updateRiskOverlay(results) {
    const probs = Object.values(results);
    const avg = probs.reduce((a, b) => a + b, 0) / probs.length;
    const pct = (avg * 100).toFixed(1);

    const el = $('#overall-risk');
    el.textContent = `${pct}%`;
    el.className = 'risk-overlay__value';
    if (avg > 0.7) el.classList.add('high');
    else if (avg > 0.4) el.classList.add('medium');
    else el.classList.add('low');

    const conf = avg > 0.6 ? 'High' : avg > 0.4 ? 'Moderate' : 'Low';
    $('#risk-confidence').textContent = `Prediction Confidence: ${conf}`;
}

/* ---------- FOOTER ---------- */
function updateFooter(smiles) {
    const now = new Date();
    const ts = now.toLocaleString('en-US', { hour: 'numeric', minute: '2-digit', second: '2-digit', hour12: true });
    $('#footer-timestamp').textContent = `Last prediction: ${ts}`;
}

function updateHistoryNav() {
    const prev = $('#history-prev');
    const next = $('#history-next');
    const label = $('#history-label');

    prev.disabled = historyIdx <= 0;
    next.disabled = historyIdx >= history.length - 1;
    label.textContent = `${historyIdx + 1} of ${history.length}`;
}

function navigateHistory(delta) {
    const newIdx = historyIdx + delta;
    if (newIdx < 0 || newIdx >= history.length) return;
    historyIdx = newIdx;
    const entry = history[historyIdx];
    smilesInput.value = entry.smiles;
    renderPredictions(entry.results);
    renderBarChart(entry.results);
    renderRadarChart(entry.results);
    renderMolecule(entry.smiles, entry.results);
    updateRiskOverlay(entry.results);
    updateHistoryNav();
}

/* ---------- CLEAR ---------- */
function handleClear() {
    smilesInput.value = '';
    emptyState.classList.remove('hidden');
    seCards.classList.add('hidden');
    chartsSection.classList.add('hidden');
    molPlaceholder.classList.remove('hidden');
    molActive.classList.add('hidden');
    riskOverlay.classList.add('hidden');
    loadingState.classList.add('hidden');
    $('#confidence-val').textContent = '—';
    $('#latency-val').textContent = '—';
    smilesInput.focus();
}

/* ---------- PASTE ---------- */
async function handlePaste() {
    try {
        const text = await navigator.clipboard.readText();
        smilesInput.value = text;
        smilesInput.focus();
    } catch {
        smilesInput.focus();
    }
}

/* ---------- UTILITY ---------- */
function getGradientColor(prob) {
    // Teal → Amber → Red gradient
    if (prob < 0.4) {
        // Teal range
        const t = prob / 0.4;
        return lerpColor('#0D9488', '#F59E0B', t);
    } else if (prob < 0.7) {
        const t = (prob - 0.4) / 0.3;
        return lerpColor('#F59E0B', '#F97316', t);
    } else {
        const t = (prob - 0.7) / 0.3;
        return lerpColor('#F97316', '#EF4444', t);
    }
}

function lerpColor(a, b, t) {
    const ah = parseInt(a.replace('#', ''), 16);
    const bh = parseInt(b.replace('#', ''), 16);
    const ar = (ah >> 16) & 0xff, ag = (ah >> 8) & 0xff, ab = ah & 0xff;
    const br = (bh >> 16) & 0xff, bg = (bh >> 8) & 0xff, bb = bh & 0xff;
    const rr = Math.round(ar + (br - ar) * t);
    const rg = Math.round(ag + (bg - ag) * t);
    const rb = Math.round(ab + (bb - ab) * t);
    return `rgb(${rr},${rg},${rb})`;
}
