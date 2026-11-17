(function() {
  var ART = 'ul[aria-labelledby="dropdown-articles"] a.dropdown-item';
  var BEN = 'ul[aria-labelledby="dropdown-benchmarks"] a.dropdown-item';
  // href$= matches both articles/<name>.html and articles/id/<name>.html
  function art(name) { return ART + '[href$="/' + name + '.html"]'; }
  function ben(name) { return BEN + '[href$="/' + name + '.html"]'; }

  var I18N = {
    nav: [
      { sel: '#navbar a.nav-link[href$="getting-started.html"]', en: 'Get Started', id: 'Mulai' },
      { sel: '#dropdown-articles', en: 'Articles', id: 'Artikel' },
      { sel: '#navbar a.nav-link[href$="reference/index.html"]', en: 'Reference', id: 'Referensi' },
      { sel: '#dropdown-benchmarks', en: 'Benchmarks', id: 'Tolok Ukur' },
      // Articles menu
      { sel: art('getting-started'), en: 'Getting Started with fastsae', id: 'Memulai fastsae' },
      { sel: art('optimal-workflow'), en: 'An Optimal SAE Workflow', id: 'Alur Kerja SAE yang Optimal' },
      { sel: art('spatial-temporal'), en: 'Spatial & Spatio-Temporal Models (EBLUP)', id: 'Model Spasial & Spatio-Temporal (EBLUP)' },
      { sel: art('hb-area-inla'), en: 'Bayesian Area-Level Models (INLA)', id: 'Model Bayesian Tingkat Area (INLA)' },
      { sel: art('hb-priors'), en: 'Default Priors of the HB Models', id: 'Prior Bawaan Model HB' },
      { sel: art('model-diagnostics'), en: 'Model Diagnostics & Calibration', id: 'Diagnostik & Kalibrasi Model' },
      { sel: art('unit-level-bhf'), en: 'Unit-Level Estimation (BHF and HB)', id: 'Estimasi Tingkat Unit (BHF dan HB)' },
      { sel: art('twofold-subarea'), en: 'Sub-Area Level Two-Fold Models', id: 'Model Subarea Dua Tingkat' },
      { sel: art('benchmarks'), en: 'Performance Benchmarks', id: 'Tolok Ukur Kinerja' },
      { sel: art('hb-area-engines'), en: 'hb_area() Engines: INLA vs JAX (CPU / Apple GPU / CUDA) vs tipsae', id: 'Engine hb_area(): INLA vs JAX (CPU / GPU Apple / CUDA) vs tipsae' },
      { sel: art('finite-population-validation'), en: 'Finite Population Validation', id: 'Validasi Populasi Terhingga' },
      { sel: art('fastsae-evaluation'), en: 'Is fastsae Better? Accuracy, Speed and Memory', id: 'Apakah fastsae Lebih Baik? Akurasi, Kecepatan, dan Memori' },
      { sel: art('glosarium'), en: 'Glossary (Bahasa Indonesia)', id: 'Daftar Istilah' },
      // Benchmarks menu
      { sel: ben('benchmarks'), en: 'fastsae vs other SAE packages', id: 'fastsae vs paket SAE lain' },
      { sel: ben('hb-area-engines'), en: 'hb_area() engines: INLA vs JAX (CPU / Apple GPU / CUDA) vs tipsae', id: 'Engine hb_area(): INLA vs JAX (CPU / GPU Apple / CUDA) vs tipsae' },
      { sel: ben('fastsae-evaluation'), en: 'Is fastsae better? Accuracy, speed and memory (SFH, spatial Beta)', id: 'Apakah fastsae lebih baik? Akurasi, kecepatan, dan memori (SFH, Beta spasial)' }
    ],
    headings: {
      "Why fastsae?": "Mengapa fastsae?",
      "Key Features": "Fitur Utama",
      "Overview": "Gambaran Umum",
      "Installation": "Instalasi",
      "Comparison with Other Packages": "Perbandingan dengan Paket Lain",
      "Supported Models": "Model yang Didukung",
      "Performance at a Glance": "Ringkasan Kinerja",
      "Documentation and User Guides": "Dokumentasi dan Panduan Pengguna",
      "References": "Daftar Pustaka"
    },
    referenceGroups: {
      "Classical Area-Level Models (EBLUP)": "Model Klasik Tingkat Area (EBLUP)",
      "Bayesian Area-Level Models (INLA)": "Model Bayesian Tingkat Area (INLA)",
      "Estimation Engines": "Engine Estimasi",
      "Two-Fold Subarea Models": "Model Subarea Dua Tingkat",
      "Unit-Level Models": "Model Tingkat Unit",
      "Pre-Modeling Utilities": "Utilitas Pra-Pemodelan",
      "Model Diagnostics & Visualization": "Diagnostik & Visualisasi Model",
      "Benchmarking & Calibration": "Benchmarking & Kalibrasi",
      "Model Comparison & Spatial Mapping": "Perbandingan Model & Pemetaan Spasial",
      "Simulation & Synthetic Data Generation": "Simulasi & Pembangkitan Data Sintetis",
      "Example & Synthetic Datasets": "Data Contoh & Data Sintetis"
    },
    referenceDescs: {
      "High-performance C++ implementation of Fay-Herriot models": "Implementasi C++ berkinerja tinggi untuk model Fay-Herriot",
      "Hierarchical Bayesian SAE with INLA (or MCMC with NumPyro / JAX) across 6 distributions with spatial and spatio-temporal effects, and sub-indicators of a parent indicator": "SAE Bayesian hierarkis dengan INLA (atau MCMC dengan NumPyro / JAX) untuk 6 distribusi, dengan efek spasial dan spatio-temporal, serta subindikator dari indikator induk",
      "Engine objects for the Hierarchical Bayes models: INLA (default) or MCMC with NumPyro / JAX (optional Python backend)": "Objek engine untuk model Bayesian hierarkis: INLA (bawaan) atau MCMC dengan NumPyro / JAX (backend Python opsional)",
      "Nested area-subarea level models (Torabi & Rao, 2014) with both area and subarea random effects via REML/ML EBLUP and hierarchical Bayes (INLA or MCMC)": "Model area-subarea tersarang (Torabi & Rao, 2014) dengan efek acak area dan subarea, melalui EBLUP REML/ML dan Bayesian hierarkis (INLA atau MCMC)",
      "Battese-Harter-Fuller models for unit-level survey data via REML/ML EBLUP and hierarchical Bayes (INLA or MCMC)": "Model Battese-Harter-Fuller untuk data survei tingkat unit, melalui EBLUP REML/ML dan Bayesian hierarkis (INLA atau MCMC)",
      "Spatial weights, spatial autocorrelation tests, and smoothing of sampling variances": "Matriks pembobot spasial, uji autokorelasi spasial, dan penghalusan varians sampling",
      "Universal diagnostic framework and ggplot2 visualization methods": "Kerangka diagnostik terpadu dan metode visualisasi ggplot2",
      "Calibrate small area estimates to aggregate targets via Ratio, Difference, Optimal, and Logit methods (Rao & Molina, 2015)": "Kalibrasi estimasi area kecil ke target agregat dengan metode rasio, selisih, optimal, dan logit (Rao & Molina, 2015)",
      "Unified spatial choropleth mapping ('satu pintu'), model concordance evaluation, and official publication exports": "Pemetaan koroplet spasial terpadu ('satu pintu'), evaluasi kesesuaian model, dan ekspor untuk publikasi resmi",
      "Utilities to simulate area-level cross-sectional, panel, and spatial weight matrices": "Utilitas untuk menyimulasikan data tingkat area lintas wilayah, data panel, dan matriks pembobot spasial",
      "Built-in empirical and synthetic survey datasets": "Data survei empiris dan sintetis bawaan"
    },
    sidebar: {
      "Links": "Tautan",
      "License": "Lisensi",
      "Community": "Komunitas",
      "Citation": "Sitasi",
      "Authors": "Penulis",
      "Developers": "Pengembang",
      "On this page": "Daftar Isi"
    },
    theme: { light: ['Light', 'Terang'], dark: ['Dark', 'Gelap'], auto: ['Auto', 'Otomatis'] }
  };

  // ---- EN <-> ID page mapping --------------------------------------------
  // Articles live in articles/<name>.html (EN) and articles/id/<name>.html (ID);
  // the home page index.html corresponds to articles/id/beranda.html, and the
  // glossary articles/id/glosarium.html has no English counterpart.
  var ID_PAGES = ['getting-started', 'optimal-workflow', 'spatial-temporal', 'hb-area-inla',
    'hb-priors', 'model-diagnostics', 'unit-level-bhf', 'twofold-subarea', 'benchmarks',
    'hb-area-engines', 'finite-population-validation', 'fastsae-evaluation'];

  function siteRoot() {
    var scripts = document.getElementsByTagName('script');
    for (var i = 0; i < scripts.length; i++) {
      var src = scripts[i].getAttribute('src') || '';
      if (/(^|\/)extra\.js(\?|$)/.test(src)) return new URL(src, document.baseURI).href.replace(/extra\.js.*$/, '');
    }
    return null;
  }

  // Path of a URL relative to the site root, or null outside the site
  function relPath(href, root) {
    var abs = new URL(href, document.baseURI).href.split('#')[0].split('?')[0];
    return abs.indexOf(root) === 0 ? abs.slice(root.length).replace(/\/$/, '') : null;
  }

  // {home|article|glossary, id} for pages with a translation, null otherwise
  function pageOf(rel) {
    if (rel === null) return null;
    if (rel === '' || rel === 'index.html') return { home: true, id: false };
    var m = rel.match(/^articles\/(id\/)?([a-z0-9-]+)\.html$/);
    if (!m) return null;
    if (m[1] && m[2] === 'beranda') return { home: true, id: true };
    if (m[1] && m[2] === 'glosarium') return { glossary: true, id: true };
    if (ID_PAGES.indexOf(m[2]) < 0) return null;
    return { article: m[2], id: !!m[1] };
  }

  function pageUrl(root, target, lang) {
    if (target.home) return root + (lang === 'id' ? 'articles/id/beranda.html' : 'index.html');
    return root + 'articles/' + (lang === 'id' ? 'id/' : '') + target.article + '.html';
  }

  // In ID mode the navigation points to the Indonesian version of each page
  function rewriteLinks(lang, root) {
    document.querySelectorAll('#navbar a[href], a.navbar-brand[href]').forEach(function(a) {
      var raw = a.getAttribute('data-orig-href');
      if (raw === null) { raw = a.getAttribute('href'); a.setAttribute('data-orig-href', raw); }
      var target = pageOf(relPath(raw, root));
      if (target && !target.glossary) {
        a.setAttribute('href', new URL(pageUrl(root, target, lang), document.baseURI).href);
      }
    });
  }

  function setText(el, text) {
    // keep icons: replace only the last text node
    for (var i = el.childNodes.length - 1; i >= 0; i--) {
      var n = el.childNodes[i];
      if (n.nodeType === 3 && n.textContent.trim() !== '') { n.textContent = ' ' + text; return; }
    }
    el.textContent = text;
  }

  function applyLanguage(lang) {
    var isId = lang === 'id';
    document.documentElement.setAttribute('data-lang', lang);
    document.documentElement.setAttribute('lang', lang);
    try { localStorage.setItem('fastsae_lang', lang); } catch (e) {}
    var root = siteRoot();
    if (root) rewriteLinks(lang, root);

    var curLabel = document.getElementById('lang-current-label');
    if (curLabel) curLabel.textContent = lang.toUpperCase();
    document.querySelectorAll('.lang-select-btn').forEach(function(btn) {
      var on = btn.getAttribute('data-lang-val') === lang;
      btn.classList.toggle('active', on);
      if (on) btn.setAttribute('aria-current', 'true'); else btn.removeAttribute('aria-current');
    });

    // 1. Navbar
    I18N.nav.forEach(function(item) {
      document.querySelectorAll(item.sel).forEach(function(el) {
        el.textContent = isId ? item.id : item.en;
      });
    });
    Object.keys(I18N.theme).forEach(function(k) {
      var b = document.querySelector('button.dropdown-item[data-bs-theme-value="' + k + '"]');
      if (b) setText(b, I18N.theme[k][isId ? 1 : 0]);
    });

    // 2. Search input
    var searchInput = document.getElementById('search-input');
    if (searchInput) {
      searchInput.setAttribute('placeholder', isId ? 'Cari fungsi / topik...' : 'Search for');
      searchInput.setAttribute('aria-label', isId ? 'Cari di situs' : 'Search site');
    }

    // 3. Headings of untranslated pages (the home page sections; the ID pages are already Indonesian)
    document.querySelectorAll('main h1, main h2, main h3').forEach(function(el) {
      if (!el.dataset.origText) el.dataset.origText = el.textContent.replace(/\s*#\s*$/, '').trim();
      var orig = el.dataset.origText;
      if (I18N.headings[orig] && el.childNodes[0]) {
        el.childNodes[0].textContent = isId ? I18N.headings[orig] : orig;
      }
    });

    // 4. Reference index groups and descriptions
    document.querySelectorAll('.section.level2 h2, .pkgdown-ref-index h2').forEach(function(h2) {
      if (!h2.dataset.origText) h2.dataset.origText = h2.textContent.trim();
      var orig = h2.dataset.origText;
      if (I18N.referenceGroups[orig]) h2.textContent = isId ? I18N.referenceGroups[orig] : orig;
    });
    document.querySelectorAll('.section.level2 p, .pkgdown-ref-index p').forEach(function(p) {
      if (!p.dataset.origText) p.dataset.origText = p.textContent.trim();
      var orig = p.dataset.origText;
      if (I18N.referenceDescs[orig]) p.textContent = isId ? I18N.referenceDescs[orig] : orig;
    });

    // 5. Sidebar and table-of-contents headings
    document.querySelectorAll('aside h2, .col-md-3 h2, nav#toc h2').forEach(function(h2) {
      if (!h2.dataset.origText) h2.dataset.origText = h2.textContent.trim();
      var orig = h2.dataset.origText;
      if (I18N.sidebar[orig]) h2.textContent = isId ? I18N.sidebar[orig] : orig;
    });

    // 6. "Source:" line of articles and the footer
    document.querySelectorAll('main small.dont-index, footer p').forEach(function(p) {
      if (!p.dataset.origHtml) p.dataset.origHtml = p.innerHTML;
      p.innerHTML = !isId ? p.dataset.origHtml : p.dataset.origHtml
        .replace(/^(\s*)Source:/, '$1Sumber:')
        .replace(/Developed by/g, 'Dikembangkan oleh')
        .replace(/Site built with/g, 'Situs dibuat dengan');
    });
  }

  // Choosing a language on a page that has a counterpart opens the counterpart
  function switchLanguage(lang) {
    var root = siteRoot();
    var cur = root ? pageOf(relPath(window.location.href, root)) : null;
    if (cur && !cur.glossary && cur.id !== (lang === 'id')) {
      try { localStorage.setItem('fastsae_lang', lang); } catch (e) {}
      window.location.href = pageUrl(root, cur, lang) + window.location.hash;
      return;
    }
    applyLanguage(lang);
  }

  function savedLang() {
    try { return localStorage.getItem('fastsae_lang') || 'en'; } catch (e) { return 'en'; }
  }

  function initLangSwitch() {
    var root = siteRoot();
    var page = root ? pageOf(relPath(window.location.href, root)) : null;
    var activeLang = savedLang();
    if (page && page.id) activeLang = 'id';   // an Indonesian page is always shown in ID mode
    // saved preference is Indonesian and this English page has a translation: go there
    if (page && !page.id && activeLang === 'id') {
      window.location.replace(pageUrl(root, page, 'id'));
      return;
    }

    if (!document.getElementById('dropdown-langswitch')) {
      var lightswitchBtn = document.getElementById('dropdown-lightswitch');
      var lightswitchLi = lightswitchBtn ? lightswitchBtn.closest('li') : null;
      var rightNav = document.querySelector('#navbar ul.navbar-nav:not(.me-auto)') ||
                     document.querySelectorAll('#navbar ul.navbar-nav')[1];
      var langLi = document.createElement('li');
      langLi.className = 'nav-item dropdown lang-switch-item';
      langLi.innerHTML =
        '<button class="nav-link dropdown-toggle" type="button" id="dropdown-langswitch" data-bs-toggle="dropdown" aria-expanded="false" aria-label="Language / Bahasa" title="Language / Bahasa">' +
          '<span class="fa fa-globe" aria-hidden="true"></span> <span id="lang-current-label">' + activeLang.toUpperCase() + '</span>' +
        '</button>' +
        '<ul class="dropdown-menu dropdown-menu-end" aria-labelledby="dropdown-langswitch">' +
          '<li><button class="dropdown-item lang-select-btn" type="button" data-lang-val="en" lang="en"><span class="me-2">🇺🇸</span> English (EN)</button></li>' +
          '<li><button class="dropdown-item lang-select-btn" type="button" data-lang-val="id" lang="id"><span class="me-2">🇮🇩</span> Bahasa Indonesia (ID)</button></li>' +
        '</ul>';

      if (lightswitchLi && lightswitchLi.parentNode) {
        lightswitchLi.parentNode.insertBefore(langLi, lightswitchLi);
      } else if (rightNav) {
        rightNav.appendChild(langLi);
      }

      langLi.querySelectorAll('.lang-select-btn').forEach(function(btn) {
        btn.addEventListener('click', function(e) {
          e.preventDefault();
          switchLanguage(this.getAttribute('data-lang-val'));
          try {
            var dropdownEl = document.getElementById('dropdown-langswitch');
            if (window.bootstrap && bootstrap.Dropdown) {
              (bootstrap.Dropdown.getInstance(dropdownEl) || new bootstrap.Dropdown(dropdownEl)).hide();
            }
          } catch (err) {}
        });
      });
    }

    applyLanguage(activeLang);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initLangSwitch);
  } else {
    initLangSwitch();
  }
})();
