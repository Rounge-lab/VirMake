(function () {
    "use strict";

    const storageKey = "crc-theme";

    function preferredTheme() {
        let savedTheme = null;

        try {
            savedTheme = localStorage.getItem(storageKey);
        } catch (error) {
            // Storage may be unavailable.
        }

        if (savedTheme === "light" || savedTheme === "dark") {
            return savedTheme;
        }

        return window.matchMedia("(prefers-color-scheme: dark)").matches
            ? "dark"
            : "light";
    }

    function applyTheme(theme) {
        document.documentElement.setAttribute("data-theme", theme);

        const button = document.getElementById("crc-theme-toggle");

        if (!button) {
            return;
        }

        const isDark = theme === "dark";
        const icon = button.querySelector("svg");

        // A stable accessible name pairs with aria-pressed:
        // pressed means dark mode is enabled.
        button.setAttribute("aria-label", "Dark mode");
        button.setAttribute("aria-pressed", String(isDark));

        button.setAttribute(
            "title",
            isDark ? "Switch to light theme" : "Switch to dark theme"
        );

        if (icon) {
            icon.setAttribute("data-theme", theme);
        }
    }

    function initializeToggle() {
        const button = document.getElementById("crc-theme-toggle");

        if (!button) {
            return;
        }

        applyTheme(
            document.documentElement.getAttribute("data-theme")
                || preferredTheme()
        );

        button.addEventListener("click", function () {
            const currentTheme =
                document.documentElement.getAttribute("data-theme");

            const nextTheme = currentTheme === "dark" ? "light" : "dark";

            try {
                localStorage.setItem(storageKey, nextTheme);
            } catch (error) {
                // The toggle still works without persistent storage.
            }

            applyTheme(nextTheme);
        });
    }

    if (document.readyState === "loading") {
        document.addEventListener("DOMContentLoaded", initializeToggle);
    } else {
        initializeToggle();
    }
}());
