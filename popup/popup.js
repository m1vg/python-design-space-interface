define([
    "nbextensions/widgets/widgets/js/widget",
    "jqueryui",
    "underscore",
    "base/js/utils",
    "bootstrap",
], function(widget, $, _, utils) {
    "use strict";

    var PopupView = widget.DOMWidgetView.extend({
        initialize: function() {
            /**
             * Public constructor
             */
            PopupView.__super__.initialize.apply(this, arguments);
            this.children_views = new widget.ViewList(this.add_child_model, null, this);
            this.listenTo(this.model, 'change:children', function(model, value) {
                this.children_views.update(value);
            }, this);
            this.listenTo(this.model, 'change:overflow_x', this.update_overflow_x, this);
            this.listenTo(this.model, 'change:overflow_y', this.update_overflow_y, this);
            this.listenTo(this.model, 'change:box_style', this.update_box_style, this);
            this._shown_once = false;
            this.popped_out = true;
        },

        update_attr: function(name, value) {
            /**
             * Set a css attr of the widget view.
             */
            this.$window.css(name, value);
        },

        render: function() {
            /**
             * Called when view is rendered.
             */
	    var that = this;
	    this.$window = $('<div style="background-color:white; border:1px solid black;" />')
                .addClass('modal widget-modal')
                .appendTo($('#notebook-container'))
                .css('border-radius', '4px')
                .css('-moz-border-radius', '4px')
                .css('-webkit-border-radius', '4px')
                .mousedown(function(){
                    that.bring_to_front();
                });
            this.$title_bar = $('<div />')
                .addClass('popover-title')
                .appendTo(this.$window)
                .mousedown(function(){
                    that.bring_to_front();
                });
            this.$title = $('<div />')
                .addClass('widget-modal-title')
                .appendTo(this.$title_bar);
                                                
            this.$close = $('<button/>')
                .css('margin-left', '3px')
                .css('-moz-box-shadow', 'inset 0px 33px 0px -24px #e67a73')
                .css('-webkit-box-shadow', 'inset 0px 33px 0px -24px #e67a73')
                .css('box-shadow', 'inset 0px 33px 0px -24px #e67a73')
                .css('background-color', '#e4685d')
                .css('border-radius', '4px')
                .css('-moz-border-radius', '4px')
                .css('-webkit-border-radius', '4px')
                .css('border', '1px solid #ffffff')
                .css('display', 'inline-block')
                .css('padding', '8px 8px')
                                                
                .appendTo(this.$title_bar)
                .click(function(){
                    that.hide();
                    event.stopPropagation();
                });
                

            this.$minimize = $('<button />')
                .css('margin-left', '5px')
                .css('-moz-box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                .css('-webkit-box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                .css('box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                .css('background-color', '#ecfc0d')
                .css('border-radius', '4px')
                .css('-moz-border-radius', '4px')
                .css('-webkit-border-radius', '4px')
                .css('border', '1px solid #ffffff')
                .css('display', 'inline-block')
                .css('padding', '8px 8px')
                .appendTo(this.$title_bar)
                .click(function(){
                    that.popped_out = !that.popped_out;
                    if (!that.popped_out) {
                        that.$minimize
                              .css('margin-left', '5px')
                              .css('-moz-box-shadow', 'inset 0px 33px 0px -24px #3dc21b')
                              .css('-webkit-box-shadow', 'inset 0px 33px 0px -24px #3dc21b')
                              .css('box-shadow', 'inset 0px 33px 0px -24px #3dc21b')
                              .css('background-color', '#44c767')
                              .css('border-radius', '4px')
                              .css('-moz-border-radius', '4px')
                              .css('-webkit-border-radius', '4px')
                              .css('border', '1px solid #ffffff')
                              .css('display', 'inline-block')
                              .css('padding', '8px 8px')
                       
                        that.$window
                            .draggable('destroy')
                            .resizable('destroy')
                            .removeClass('widget-modal modal')
                            .addClass('docked-widget-modal')
                            .detach()

                        that.$title_bar
                            .detach()
                            .appendTo(that.$el)
                        that.$body
                            .detach()
                            .appendTo(that.$el)
                        that.$show_button.hide();
                        that.$close.hide();
                    } else {
                        that.$minimize
                               .css('margin-left', '5px')
                               .css('-moz-box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                               .css('-webkit-box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                               .css('box-shadow', 'inset 0px 33px 0px -24px #ffff00')
                               .css('background-color', '#ecfc0d')
                               .css('border-radius', '4px')
                               .css('-moz-border-radius', '4px')
                               .css('-webkit-border-radius', '4px')
                               .css('border', '1px solid #ffffff')
                               .css('display', 'inline-block')
                               .css('padding', '8px 8px')

                        that.$window
                            .css('border', '1px solid black')
                            .removeClass('docked-widget-modal')
                            .addClass('widget-modal modal')
                            .appendTo($('#notebook-container'))
                            .draggable({handle: '.popover-title',
                                       snap: '#notebook, .modal',
                                       snapMode: 'both',
                                       })
                            .resizable()
                            .children('.ui-resizable-handle').show();
                        that.$title_bar
                            .detach()
                            .appendTo(that.$window)
                        that.$body
                            .detach()
                            .appendTo(that.$window)
                        that.show();
                        that.$show_button.show();
                        that.$close.show();
                    }
                    event.stopPropagation();
                });
            
            this.$show_button = $('<button />')
                .addClass('btn btn-info widget-modal-show')
                .appendTo(this.$el)
                .click(function(){
                    that.show();
                });
            
            this.$window.draggable({handle: '.popover-title', snap: '#notebook, .modal', snapMode: 'both'});
            this.$window.resizable();
            this.$window.on('resize', function(){
                that.$window
                    .css('position', 'fixed')
                    .css('margin-left', that.$window.off.left)
                    .css('margin-top', that.$window.off.top)
                    .css('border', '1px solid black')
                that.$body.outerHeight(that.$window.innerHeight() - that.$title_bar.outerHeight());
                that.$body.outerWidth(that.$window.innerWidth());
            });
            this.$window.on('start', function(){
                that.$window.off = that.$window.offset();
            });
            
            this.$body = $('<div />')
                .addClass('modal-body')
                .addClass('widget-modal-body')
                .addClass('widget-container')
                .addClass('vbox')

                                                
            this.$el_to_style = this.$body;
            this.$box = this.$body;
            this.$box.addClass('widget-box');
            this.children_views.update(this.model.get('children'));
            this.$body
                .appendTo(this.$window)
                .css('overflow-x', 'auto')
                .css('overflow-y', 'auto')
            
            this.update();
            this.update_overflow_x();
            this.update_overflow_y();
            this.update_box_style('');
        },

	hide: function() {
            // Called when the modal hide button is clicked.
            this.$window.hide();
            this.$show_button.removeClass('btn-info');
        },
        
    show: function() {
            // Called when the modal show button is clicked.
            this.$show_button.addClass('btn-info');
            this.$window.show();
            this.$body.show();
            if (this.popped_out) {
                this.$window
                    .css("positon", "fixed")
                    .css("width", "45%")
                    .css("height", "45%")
                    .css("top", "200px")
                    .css('border', '1px solid black')
                
                this.$window.css("left", Math.max(0, (($('body').outerWidth() - this.$window.outerWidth()) / 2) +
                    $(window).scrollLeft()) + "px");
                this.$body.outerHeight(this.$window.innerHeight() - this.$title_bar.outerHeight());
                this.$body.outerWidth(this.$window.innerWidth());
                this.bring_to_front();
            }
        },

	bring_to_front: function() {
            // Make the modal top-most, z-ordered about the other modals.
            var $widget_modals = $(".widget-modal");
            var max_zindex = 0;
            $widget_modals.each(function (index, el){
                max_zindex = Math.max(max_zindex, parseInt($(el).css('z-index')));
            });
            
            // Start z-index of widget modals at 2000
            max_zindex = Math.max(max_zindex, 2000);

            $widget_modals.each(function (index, el){
                var $el = $(el);
                if (max_zindex == parseInt($el.css('z-index'))) {
                    $el.css('z-index', max_zindex - 1);
                }
            });
            this.$window.css('z-index', max_zindex);
        },

        update: function(){
                if (!this._shown_once) {
                        this._shown_once = true;
                        this.show();
                }
                return PopupView.__super__.update.apply(this);
        },
        update_overflow_x: function() {
            /**
             * Called when the x-axis overflow setting is changed.
             */
            this.$box.css('overflow-x', this.model.get('overflow_x'));
        },

        update_overflow_y: function() {
            /**
             * Called when the y-axis overflow setting is changed.
             */
            this.$box.css('overflow-y', this.model.get('overflow_y'));
        },
        
        update_box_style: function(previous_trait_value) {
            var class_map = {
                success: ['alert', 'alert-success'],
                info: ['alert', 'alert-info'],
                warning: ['alert', 'alert-warning'],
                danger: ['alert', 'alert-danger']
            };
            this.update_mapped_classes(class_map, 'box_style', previous_trait_value, this.$box[0]);
        },

        add_child_model: function(model) {
            /**
             * Called when a model is added to the children list.
             */
            var that = this;
            var dummy = $('<div/>');
            that.$box.append(dummy);
            var view =  this.create_child_view(model).then(function(view) {
                dummy.replaceWith(view.el);
                // Trigger the displayed event of the child view.
                that.displayed.then(function() {
                    view.trigger('displayed');
                });
                return view;
            }).catch(utils.reject("Couldn't add child view to box", true));
            return view;
        },

        remove: function() {
            /**
             * We remove this widget before removing the children as an optimization
             * we want to remove the entire container from the DOM first before
             * removing each individual child separately.
             */
            PopupView.__super__.remove.apply(this, arguments);
            this.children_views.remove();
            this.$window.remove();
        },
        
        _get_selector_element: function(selector) {
            // Get an element view a 'special' jquery selector.  (see widget.js)
            //
            // Since the modal actually isn't within the $el in the DOM, we need to extend
            // the selector logic to allow the user to set css on the modal if need be.
            // The convention used is:
            // "modal" - select the modal div
            // "modal [selector]" - select element(s) within the modal div.
            // "[selector]" - select elements within $el
            // "" - select the $el_to_style
            if (selector.substring(0, 5) == 'modal') {
                if (selector == 'modal') {
                    return this.$window;
                } else {
                    return this.$window.find(selector.substring(6));
                }
            } else {
                return PopupView.__super__._get_selector_element.apply(this, [selector]);
            }
        },
        
    });



    return {
        PopupView: PopupView
    }
});