/* ================================================================
   MedAI Insights – Landing Page Interactions
   ================================================================ */

document.addEventListener('DOMContentLoaded', () => {
  initNavbar();
  initParticles();
  initMiniChart();
  initScrollReveal();
  initCountUp();
});

/* ---------- Navbar scroll effect ---------- */
function initNavbar() {
  const navbar = document.getElementById('navbar');
  const toggle = document.getElementById('nav-toggle');
  const links = document.getElementById('nav-links');

  window.addEventListener('scroll', () => {
    navbar.classList.toggle('scrolled', window.scrollY > 40);
  }, { passive: true });

  toggle.addEventListener('click', () => {
    links.classList.toggle('active');
  });

  // Close mobile menu on link click
  links.querySelectorAll('a').forEach(a => {
    a.addEventListener('click', () => links.classList.remove('active'));
  });
}

/* ---------- Particle canvas ---------- */
function initParticles() {
  const canvas = document.getElementById('particles');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  function resize() {
    canvas.width = canvas.parentElement.offsetWidth;
    canvas.height = canvas.parentElement.offsetHeight;
  }
  resize();
  window.addEventListener('resize', resize);

  const PARTICLE_COUNT = 60;
  const particles = [];
  const colors = ['rgba(13,148,136,0.30)', 'rgba(6,182,212,0.25)', 'rgba(20,184,166,0.22)', 'rgba(16,185,129,0.18)'];

  for (let i = 0; i < PARTICLE_COUNT; i++) {
    particles.push({
      x: Math.random() * canvas.width,
      y: Math.random() * canvas.height,
      r: Math.random() * 2.2 + 0.6,
      dx: (Math.random() - 0.5) * 0.3,
      dy: (Math.random() - 0.5) * 0.3,
      color: colors[Math.floor(Math.random() * colors.length)],
    });
  }

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    for (const p of particles) {
      p.x += p.dx;
      p.y += p.dy;
      if (p.x < 0) p.x = canvas.width;
      if (p.x > canvas.width) p.x = 0;
      if (p.y < 0) p.y = canvas.height;
      if (p.y > canvas.height) p.y = 0;

      ctx.beginPath();
      ctx.arc(p.x, p.y, p.r, 0, Math.PI * 2);
      ctx.fillStyle = p.color;
      ctx.fill();
    }

    // Connection lines
    for (let i = 0; i < particles.length; i++) {
      for (let j = i + 1; j < particles.length; j++) {
        const dx = particles[i].x - particles[j].x;
        const dy = particles[i].y - particles[j].y;
        const dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < 140) {
          ctx.beginPath();
          ctx.moveTo(particles[i].x, particles[i].y);
          ctx.lineTo(particles[j].x, particles[j].y);
          ctx.strokeStyle = `rgba(13,148,136,${0.08 * (1 - dist / 140)})`;
          ctx.lineWidth = 0.7;
          ctx.stroke();
        }
      }
    }
    requestAnimationFrame(draw);
  }
  draw();
}

/* ---------- Mini probability chart bars ---------- */
function initMiniChart() {
  const chart = document.getElementById('mini-chart');
  if (!chart) return;
  const heights = [55, 75, 45, 90, 60, 82, 50, 70, 88, 40];
  heights.forEach(h => {
    const bar = document.createElement('span');
    bar.style.height = (h * 0.45) + 'px';
    chart.appendChild(bar);
  });
}

/* ---------- Scroll reveal (Intersection Observer) ---------- */
function initScrollReveal() {
  const reveals = document.querySelectorAll('.reveal');
  if (!reveals.length) return;

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        entry.target.classList.add('visible');
        observer.unobserve(entry.target);
      }
    });
  }, { threshold: 0.15 });

  reveals.forEach(el => observer.observe(el));
}

/* ---------- Count-Up Animation ---------- */
function initCountUp() {
  const numbers = document.querySelectorAll('.stat-card__number[data-count]');
  if (!numbers.length) return;

  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        animateCount(entry.target);
        observer.unobserve(entry.target);
      }
    });
  }, { threshold: 0.5 });

  numbers.forEach(el => observer.observe(el));
}

function animateCount(el) {
  const target = parseFloat(el.dataset.count);
  const suffix = el.dataset.suffix || '';
  const isDecimal = el.dataset.decimal === 'true';
  const duration = 2000;
  const start = performance.now();

  function tick(now) {
    const progress = Math.min((now - start) / duration, 1);
    const eased = 1 - Math.pow(1 - progress, 3); // ease-out cubic
    let current = eased * target;

    if (isDecimal) {
      el.textContent = current.toFixed(1) + suffix;
    } else if (target >= 1000000) {
      el.textContent = (current / 1000000).toFixed(1) + 'M' + suffix;
    } else {
      el.textContent = Math.floor(current).toLocaleString() + suffix;
    }

    if (progress < 1) requestAnimationFrame(tick);
  }
  requestAnimationFrame(tick);
}



