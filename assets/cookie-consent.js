// Cookie Consent Manager for Matomo Tracking
(function() {
  'use strict';
  
  const COOKIE_NAME = 'matomo_consent';
  const COOKIE_EXPIRY_DAYS = 365;
  
  // Get cookie by name
  function getCookie(name) {
    const value = `; ${document.cookie}`;
    const parts = value.split(`; ${name}=`);
    if (parts.length === 2) return parts.pop().split(';').shift();
    return null;
  }
  
  // Set cookie
  function setCookie(name, value, days) {
    const date = new Date();
    date.setTime(date.getTime() + (days * 24 * 60 * 60 * 1000));
    const expires = `expires=${date.toUTCString()}`;
    document.cookie = `${name}=${value};${expires};path=/;SameSite=Lax`;
  }
  
  // Load Matomo tracking script
  function loadMatomo() {
    if (window._paq) {
      console.log('Matomo already loaded');
      return;
    }
    
    window._paq = window._paq || [];
    _paq.push(['trackPageView']);
    _paq.push(['enableLinkTracking']);
    
    const u = "https://matomo.sib.swiss/";
    _paq.push(['setTrackerUrl', u + 'matomo.php']);
    _paq.push(['setSiteId', '220']);
    
    const d = document;
    const g = d.createElement('script');
    const s = d.getElementsByTagName('script')[0];
    g.async = true;
    g.src = u + 'matomo.js';
    s.parentNode.insertBefore(g, s);
    
    console.log('Matomo tracking loaded');
  }
  
  // Create and show cookie banner
  function showCookieBanner() {
    const banner = document.createElement('div');
    banner.id = 'cookie-consent-banner';
    banner.innerHTML = `
      <div class="cookie-consent-content">
        <p class="cookie-consent-text">
          We use cookies to analyze website traffic and improve your experience. 
          By accepting, you allow us to use analytics cookies.
        </p>
        <div class="cookie-consent-buttons">
          <button id="cookie-accept" class="cookie-btn cookie-btn-accept">Accept</button>
          <button id="cookie-decline" class="cookie-btn cookie-btn-decline">Decline</button>
        </div>
      </div>
    `;
    
    document.body.appendChild(banner);
    
    // Handle accept button
    document.getElementById('cookie-accept').addEventListener('click', function() {
      setCookie(COOKIE_NAME, 'accepted', COOKIE_EXPIRY_DAYS);
      banner.remove();
      loadMatomo();
    });
    
    // Handle decline button
    document.getElementById('cookie-decline').addEventListener('click', function() {
      setCookie(COOKIE_NAME, 'declined', COOKIE_EXPIRY_DAYS);
      banner.remove();
    });
  }
  
  // Initialize on page load
  function init() {
    const consent = getCookie(COOKIE_NAME);
    
    if (consent === 'accepted') {
      loadMatomo();
    } else if (consent !== 'declined') {
      // Show banner if no decision has been made
      showCookieBanner();
    }
  }
  
  // Run when DOM is ready
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
