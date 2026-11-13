(function () {
    "use strict";

    const storageKey = "crc-theme";

    function preferredTheme() {
        const savedTheme = localStorage.getItem(storageKey);

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

        if (button) {
            const isDark = theme === "dark";

            button.setAttribute(
                "aria-label",
                isDark ? "Switch to light theme" : "Switch to dark theme"
            );

            button.setAttribute(
                "title",
                isDark ? "Switch to light theme" : "Switch to dark theme"
            );

            button.textContent = isDark ? "☀" : "☾";
        }
    }

    function addToggleButton() {
        const button = document.createElement("button");

        button.id = "crc-theme-toggle";
        button.type = "button";
        button.className = "crc-theme-toggle";

        button.addEventListener("click", function () {
            const currentTheme =
                document.documentElement.getAttribute("data-theme");

            const nextTheme = currentTheme === "dark" ? "light" : "dark";

            localStorage.setItem(storageKey, nextTheme);
            applyTheme(nextTheme);
        });

        document.body.appendChild(button);
        applyTheme(document.documentElement.getAttribute("data-theme"));
    }

    document.addEventListener("DOMContentLoaded", function () {
        applyTheme(preferredTheme());
        addToggleButton();
    });
}());
