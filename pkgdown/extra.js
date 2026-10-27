(function() {
  var I18N = {
    nav: [
      { sel: '#navbar a[href*="getting-started.html"]', en: 'Get Started', id: 'Mulai' },
      { sel: '#dropdown-articles', en: 'Articles', id: 'Artikel' },
      { sel: '#navbar a[href*="reference/index.html"], #navbar a[href*="reference/"]', en: 'Reference', id: 'Referensi' },
      { sel: '#navbar a[href*="benchmarks.html"]', en: 'Benchmarks', id: 'Benchmark' },
      // Dropdown articles (technical terms stay in English, as in the Indonesian articles)
      { sel: 'a.dropdown-item[href*="getting-started.html"]', en: 'Getting Started with fastsae', id: 'Memulai dengan fastsae' },
      { sel: 'a.dropdown-item[href*="spatial-temporal.html"]', en: 'Spatial & Spatio-Temporal Models (EBLUP)', id: 'Model Spatial & Spatio-Temporal (EBLUP)' },
      { sel: 'a.dropdown-item[href*="hb-area-inla.html"]', en: 'Bayesian Area-Level Models (INLA)', id: 'Model Bayesian Area-Level (INLA)' },
      { sel: 'a.dropdown-item[href*="hb-priors.html"]', en: 'Default Priors of the HB Models', id: 'Default Prior Model HB' },
      { sel: 'a.dropdown-item[href*="model-diagnostics.html"]', en: 'Model Diagnostics & Calibration', id: 'Diagnostik Model & Kalibrasi' },
      { sel: 'a.dropdown-item[href*="unit-level-bhf.html"]', en: 'Unit-Level Estimation (BHF and HB)', id: 'Estimasi Unit-Level (BHF dan HB)' },
      { sel: 'a.dropdown-item[href*="twofold-subarea.html"]', en: 'Sub-Area Level Two-Fold Models', id: 'Model Two-Fold Sub-Area Level' },
      { sel: '#navbar .dropdown-menu a.dropdown-item[href$="articles/benchmarks.html"], #navbar .dropdown-menu a.dropdown-item[href$="articles/id/benchmarks.html"]', en: 'Performance Benchmarks', id: 'Benchmark Kinerja' },
      { sel: 'a.dropdown-item[href*="hb-area-engines.html"]', en: 'hb_area() Engines: INLA vs JAX (CPU / Apple GPU / CUDA) vs tipsae', id: 'Engine hb_area(): INLA vs JAX (CPU / Apple GPU / CUDA) vs tipsae' },
      { sel: 'a.dropdown-item[href*="finite-population-validation.html"]', en: 'Finite Population Validation', id: 'Validasi Finite Population' }
    ],
    headings: {
      "Why fastsae?": "Mengapa fastsae?",
      "Key Features": "Fitur Utama",
      "Overview": "Gambaran Umum",
      "Installation": "Instalasi",
      "Exact Numerical Equivalence with sae": "Ekuivalensi Numerik Eksak dengan sae",
      "Comparison with Other Packages": "Perbandingan dengan Paket Lain",
      "Model Diagnostic System": "Sistem Diagnostik Model",
      "Hierarchical Two-Stage Benchmarking": "Benchmarking Hierarkis Dua Tahap",
      "Model Concordance Evaluation": "Evaluasi Keselarasan Model",
      "Spatial Mapping API": "API Pemetaan Spasial (\"Satu Pintu\")",
      "Quick Start: Classical Fay-Herriot (EBLUP)": "Panduan Cepat: Fay-Herriot Klasik (EBLUP)",
      "Quick Start: Bayesian Area-Level (INLA)": "Panduan Cepat: Bayesian Area-Level (INLA)",
      "Performance Benchmarks": "Tolok Ukur Kinerja Komputasi",
      "Finite Population Validation": "Validasi Simulasi Populasi Terhingga",
      "Getting Started with fastsae": "Memulai Penggunaan fastsae",
      "Unit-Level Estimation (BHF)": "Estimasi Unit-Level (BHF)",
      "Spatial & Spatio-Temporal Models (EBLUP)": "Model Spasial & Spatio-Temporal (EBLUP)",
      "Bayesian Area-Level Models (INLA)": "Model Bayesian Area-Level (INLA)",
      "Model Diagnostics & Calibration": "Diagnostik Model & Kalibrasi"
    },
    referenceGroups: {
      "Classical Area-Level Models (EBLUP)": "Model Klasik Area-Level (EBLUP)",
      "Bayesian Area-Level Models (INLA / HB)": "Model Bayesian Area-Level (INLA / HB)",
      "Unit-Level Models": "Model Unit-Level",
      "Model Diagnostics & Visualization": "Diagnostik Model & Visualisasi",
      "Benchmarking & Calibration": "Benchmarking & Kalibrasi",
      "Model Comparison & Spatial Mapping": "Perbandingan Model & Pemetaan Spasial",
      "Simulation & Synthetic Data Generation": "Simulasi & Pembangkitan Data Sintetis",
      "Example & Synthetic Datasets": "Dataset Contoh & Sintetis Bawaan"
    },
    referenceDescs: {
      "High-performance C++ implementation of Fay-Herriot models": "Implementasi C++ berkinerja tinggi untuk model Fay-Herriot",
      "Hierarchical Bayesian Small Area Estimation using INLA across 6 distributions with spatial and spatio-temporal effects": "Estimasi Area Kecil Hierarchical Bayesian menggunakan INLA pada 6 distribusi dengan efek spasial dan spatio-temporal",
      "Battese-Harter-Fuller model for unit-level survey data": "Model Battese-Harter-Fuller untuk data survei tingkat unit/individu",
      "Universal diagnostic framework and ggplot2 visualization methods": "Kerangka diagnostik universal dan metode visualisasi ggplot2",
      "Calibrate small area estimates to aggregate targets via Ratio, Difference, Optimal, and Logit methods (Rao & Molina, 2015)": "Kalibrasi estimasi area kecil ke target agregat via metode Rasio, Selisih, Optimal, dan Logit (Rao & Molina, 2015)",
      "Unified spatial choropleth mapping ('satu pintu'), model concordance evaluation, and official publication exports": "Pemetaan choropleth spasial terpadu ('satu pintu'), evaluasi keselarasan model, dan ekspor publikasi resmi",
      "Utilities to simulate area-level cross-sectional, panel, and spatial weight matrices": "Utilitas untuk simulasi data cross-section area, panel, dan matriks bobot spasial",
      "Built-in empirical and synthetic survey datasets": "Dataset survei empiris dan sintetis bawaan"
    },
    sidebar: {
      "Links": "Tautan",
      "License": "Lisensi",
      "Community": "Komunitas",
      "Citation": "Sitasi",
      "Authors": "Penulis",
      "Developers": "Pengembang"
    }
  };

  // ---- EN <-> ID page mapping --------------------------------------------
  // Articles live in articles/<name>.html (EN) and articles/id/<name>.html (ID);
  // the home page index.html corresponds to articles/id/beranda.html.
  var ID_PAGES = ['getting-started', 'spatial-temporal', 'hb-area-inla', 'hb-priors',
    'model-diagnostics', 'unit-level-bhf', 'twofold-subarea', 'benchmarks',
    'hb-area-engines', 'finite-population-validation'];

  function siteRoot() {
    var scripts = document.getElementsByTagName('script');
    for (var i = 0; i < scripts.length; i++) {
      var src = scripts[i].getAttribute('src') || '';
      if (/(^|\/)extra\.js(\?|$)/.test(src)) return new URL(src, document.baseURI).href.replace(/extra\.js.*$/, '');
    }
    return null;
  }

  // 'home', {article: name, id: bool} or null for pages that have no translation
  function currentPage() {
    var root = siteRoot();
    if (!root) return null;
    var path = window.location.href.split('#')[0].split('?')[0];
    if (path.indexOf(root) !== 0) return null;
    var rel = path.slice(root.length).replace(/\/$/, '');
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

  // In ID mode the navigation points to the Indonesian version of each article
  function rewriteLinks(lang) {
    var root = siteRoot();
    if (!root) return;
    document.querySelectorAll('#navbar a[href], .navbar-brand[href]').forEach(function(a) {
      var raw = a.getAttribute('data-orig-href');
      if (raw === null) { raw = a.getAttribute('href'); a.setAttribute('data-orig-href', raw); }
      var abs = new URL(raw, document.baseURI).href.split('#')[0];
      var rel = abs.indexOf(root) === 0 ? abs.slice(root.length) : null;
      if (rel === null) return;
      var target = null;
      if (rel === '' || rel === 'index.html') target = { home: true };
      else {
        var m = rel.match(/^articles\/([a-z0-9-]+)\.html$/);
        if (m && ID_PAGES.indexOf(m[1]) >= 0) target = { article: m[1] };
      }
      if (target) a.setAttribute('href', new URL(pageUrl(root, target, lang), document.baseURI).href);
    });
  }

  function switchLanguage(lang) {
    var cur = currentPage();
    var root = siteRoot();
    if (cur && root && !cur.glossary && cur.id !== (lang === 'id')) {
      localStorage.setItem('fastsae_lang', lang);
      window.location.href = pageUrl(root, cur, lang);
      return;
    }
    applyLanguage(lang);
  }

  function applyLanguage(lang) {
    document.documentElement.setAttribute('data-lang', lang);
    document.documentElement.setAttribute('lang', lang);
    localStorage.setItem('fastsae_lang', lang);
    rewriteLinks(lang);

    var curLabel = document.getElementById('lang-current-label');
    if (curLabel) curLabel.textContent = lang.toUpperCase();

    document.querySelectorAll('.lang-select-btn').forEach(function(btn) {
      if (btn.getAttribute('data-lang-val') === lang) {
        btn.classList.add('active');
        btn.setAttribute('aria-current', 'true');
      } else {
        btn.classList.remove('active');
        btn.removeAttribute('aria-current');
      }
    });

    // 1. Translate nav elements
    I18N.nav.forEach(function(item) {
      document.querySelectorAll(item.sel).forEach(function(el) {
        if (!el.dataset.origText) el.dataset.origText = el.textContent.trim();
        el.textContent = (lang === 'id') ? item.id : item.en;
      });
    });

    // 2. Search input
    var searchInput = document.getElementById('search-input');
    if (searchInput) {
      if (lang === 'id') {
        searchInput.setAttribute('placeholder', 'Cari fungsi / topik...');
        searchInput.setAttribute('aria-label', 'Cari di situs');
      } else {
        searchInput.setAttribute('placeholder', 'Search for');
        searchInput.setAttribute('aria-label', 'Search site');
      }
    }

    // 3. TOC sidebar
    document.querySelectorAll('nav#toc h2, .toc-section h2, aside nav h2').forEach(function(h2) {
      if (!h2.dataset.origText) h2.dataset.origText = h2.textContent;
      h2.textContent = (lang === 'id') ? 'Daftar Isi' : h2.dataset.origText;
    });

    // 4. Headings
    document.querySelectorAll('main h1, main h2, main h3').forEach(function(el) {
      var raw = el.textContent.trim();
      var clean = raw.replace(/\s*#\s*$/, '').trim();
      if (!el.dataset.origText) el.dataset.origText = clean;
      var orig = el.dataset.origText;
      if (I18N.headings[orig]) {
        if (lang === 'id') {
          el.childNodes[0].textContent = I18N.headings[orig];
        } else {
          el.childNodes[0].textContent = orig;
        }
      }
    });

    // 5. Reference groups & descriptions
    document.querySelectorAll('.section.level2 h2, .pkgdown-ref-index h2').forEach(function(h2) {
      var t = h2.textContent.trim();
      if (!h2.dataset.origText) h2.dataset.origText = t;
      var orig = h2.dataset.origText;
      if (I18N.referenceGroups[orig]) {
        h2.textContent = (lang === 'id') ? I18N.referenceGroups[orig] : orig;
      }
    });
    document.querySelectorAll('.section.level2 p, .pkgdown-ref-index p').forEach(function(p) {
      var t = p.textContent.trim();
      if (!p.dataset.origText) p.dataset.origText = t;
      var orig = p.dataset.origText;
      if (I18N.referenceDescs[orig]) {
        p.textContent = (lang === 'id') ? I18N.referenceDescs[orig] : orig;
      }
    });

    // 6. Sidebar headings
    document.querySelectorAll('aside h2, .col-md-3 h2').forEach(function(h2) {
      var t = h2.textContent.trim();
      if (!h2.dataset.origText) h2.dataset.origText = t;
      var orig = h2.dataset.origText;
      if (I18N.sidebar[orig]) {
        h2.textContent = (lang === 'id') ? I18N.sidebar[orig] : orig;
      }
    });

    // 7. Footer
    document.querySelectorAll('footer p').forEach(function(p) {
      if (!p.dataset.origHtml) p.dataset.origHtml = p.innerHTML;
      if (lang === 'id') {
        p.innerHTML = p.dataset.origHtml
          .replace(/Developed by/g, 'Dikembangkan oleh')
          .replace(/Site built with/g, 'Situs dibuat menggunakan');
      } else {
        p.innerHTML = p.dataset.origHtml;
      }
    });
  }

  function initLangSwitch() {
    var lightswitchBtn = document.getElementById('dropdown-lightswitch');
    var lightswitchLi = lightswitchBtn ? lightswitchBtn.closest('li') : null;
    var rightNav = document.querySelector('#navbar ul.navbar-nav:not(.me-auto)') ||
                   document.querySelector('#navbar ul.navbar-nav.ms-auto') ||
                   document.querySelectorAll('#navbar ul.navbar-nav')[1];

    if (!document.getElementById('dropdown-langswitch')) {
      var curLang = localStorage.getItem('fastsae_lang') || 'en';
      var langLi = document.createElement('li');
      langLi.className = 'nav-item dropdown lang-switch-item';
      langLi.innerHTML =
        '<button class="nav-link dropdown-toggle" type="button" id="dropdown-langswitch" data-bs-toggle="dropdown" aria-expanded="false" aria-label="Language switch" title="Switch Language">' +
          '<span class="fa fa-globe"></span> <span id="lang-current-label">' + curLang.toUpperCase() + '</span>' +
        '</button>' +
        '<ul class="dropdown-menu dropdown-menu-end" aria-labelledby="dropdown-langswitch">' +
          '<li><button class="dropdown-item lang-select-btn' + (curLang === 'en' ? ' active' : '') + '" type="button" data-lang-val="en"><span class="me-2">🇺🇸</span> English (EN)</button></li>' +
          '<li><button class="dropdown-item lang-select-btn' + (curLang === 'id' ? ' active' : '') + '" type="button" data-lang-val="id"><span class="me-2">🇮🇩</span> Bahasa Indonesia (ID)</button></li>' +
        '</ul>';

      if (lightswitchLi && lightswitchLi.parentNode) {
        lightswitchLi.parentNode.insertBefore(langLi, lightswitchLi);
      } else if (rightNav) {
        rightNav.appendChild(langLi);
      }

      langLi.querySelectorAll('.lang-select-btn').forEach(function(btn) {
        btn.addEventListener('click', function(e) {
          e.preventDefault();
          var chosen = this.getAttribute('data-lang-val');
          switchLanguage(chosen);
          try {
            var dropdownEl = document.getElementById('dropdown-langswitch');
            if (window.bootstrap && bootstrap.Dropdown) {
              var inst = bootstrap.Dropdown.getInstance(dropdownEl) || new bootstrap.Dropdown(dropdownEl);
              inst.hide();
            }
          } catch(err) {}
        });
      });
    }

    var activeLang = localStorage.getItem('fastsae_lang') || 'en';
    var page = currentPage();
    if (page && page.id) activeLang = 'id';   // an Indonesian article is always shown in ID mode
    // saved preference is Indonesian but this page has a translation: go there
    if (page && !page.id && !page.glossary && activeLang === 'id' && siteRoot()) {
      window.location.replace(pageUrl(siteRoot(), page, 'id'));
      return;
    }
    applyLanguage(activeLang);
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initLangSwitch);
  } else {
    initLangSwitch();
  }
})();
